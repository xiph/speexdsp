#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "speex/speex_resampler.h"

#include <unistd.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <vector>

/**
 *  \brief A simple cli tool that can be used to read audio from a pcm file (the
 *  parameters are generally assumed to be 1ch @ 44.1kHz in int16_t encoding),
 *  resample the audio stream to another frequency (typically 16kHz).  An output
 *  pcm file will be generated.  This test harness can be used to measure the
 *  performance of the resampler with arbitrary inputs.
 *
 *  \param arg[0] -> program name
 *  \param arg[1] -> input.pcm
 *  \param arg[2] -> output.pcm
 *  \param arg[3] ... options
 *
 *  Options:
 *      -i: input sampling rate (default is 44100Hz)
 *      -o: output sampling rate (default is 16000)
 *      -c: number of channels (default is 1)
 *      -n: frame size in samples (default is 1764 for 40ms frames at 44.1kHz
 *      -q: resampler quality. Higher is 'better' with a range of [0, 10]
 * (inclusive) -h: display help message
 */

static constexpr uint32_t kMaxSamplingRate = 192000;
static constexpr uint32_t kMaxChannels = 8;
static constexpr uint32_t kMaxFrameSize = 8192;

void usage() {
  std::cout << "REQUIRED: arg[0] -> program name \n";
  std::cout << "REQUIRED: arg[1] -> input.pcm \n";
  std::cout << "REQUIRED: arg[2] -> output.pcm \n";
  std::cout << "Options: \n";
  std::cout << "  -i: input sampling rate (default is 44100Hz)\n";
  std::cout << "  -o: output sampling rate (default is 16000)\n";
  std::cout << "  -c: number of channels (default is 1)\n";
  std::cout << "  -n: frame size in samples (default is 1764 for 40ms frames "
               "at 44.1kHz\n";
  std::cout << "  -q: resampler quality. Default = 4. Higher is 'better' with "
               "a range of [0, 10] (inclusive)\n";
  std::cout << "  -h: display this message\n";
}

void print_resampler_params(int num_channels, int input_rate, int output_rate,
                            int quality) {
  printf("------ Resampler Parameters ------\n");
  printf(
      "\tChannels = %d\n\tInput Rate = %d\n\tOutput Rate = %d\n\tQuality = "
      "%d\n",
      num_channels, input_rate, output_rate, quality);
  printf("----------------------------------\n");
}

int main(int argc, char *const argv[]) {
  /* Default Options */
  uint32_t input_rate = 44100, output_rate = 16000, num_channels = 1;
  uint32_t input_frame_samples = 1764,
           resampler_quality = SPEEX_RESAMPLER_QUALITY_DEFAULT;
  std::string input_file_name, output_file_name;

  if (argc < 3) {
    std::cerr << "ERROR: Not enough input arguments";
    usage();
    return -1;
  }

  int c = 0;
  while ((c = getopt(argc, argv, "i:o:c:n:q:h:")) != -1) {
    switch (c) {
      case 'i':
        input_rate = atoi(optarg);
        break;
      case 'o':
        output_rate = atoi(optarg);
        break;
      case 'c':
        num_channels = atoi(optarg);
        break;
      case 'n':
        input_frame_samples = atoi(optarg);
        break;
      case 'q':
        resampler_quality = atoi(optarg);
        break;
      case 'h':
        usage();
        return -1;
      default:
        std::cerr << "Unidentified option: " << c;
        usage();
        return -1;
    }
  }

  // optind is for the extra arguments which are not parsed. First should point
  // to the input filename and the second to the output file name
  for (int m = 0, n = optind; n < argc && m < 2; ++m, ++n) {
    if (m == 0) {
      std::copy(argv[n], argv[n] + strlen(argv[n]),
                std::back_inserter(input_file_name));
    } else {
      std::copy(argv[n], argv[n] + strlen(argv[n]),
                std::back_inserter(output_file_name));
    }
  }

  if (input_rate > kMaxSamplingRate) {
    std::cerr << "Erroneous input rate.  Input Rate must be 0 < fs < "
              << kMaxSamplingRate;
    return -1;
  }

  if (output_rate > kMaxSamplingRate) {
    std::cerr << "Erroneous output rate.  output Rate must be 0 < fs < "
              << kMaxSamplingRate;
    return -1;
  }

  if (num_channels > kMaxChannels) {
    std::cerr << "Erroneous number of channels. Must be 0 < ch < "
              << kMaxChannels;
    return -1;
  }

  if (input_frame_samples > kMaxFrameSize) {
    std::cerr << "Erroneous frame size. Must be 0 < N < " << kMaxFrameSize;
    return -1;
  }

  if (resampler_quality > SPEEX_RESAMPLER_QUALITY_MAX) {
    std::cerr << "Erroneous output rate.  output Rate must be 0 < fs < "
              << kMaxSamplingRate;
    return -1;
  }

  std::ifstream input_file(input_file_name, std::ios::binary);
  if (!input_file.is_open()) {
    std::cerr << "Unable to open input file " << argv[1] << "\n";
    return -1;
  }
  std::cout << "Opening input file = " << input_file_name << "\n";

  std::ofstream output_file(output_file_name, std::ios::binary);
  if (!output_file.is_open()) {
    std::cerr << "Unable to open output file " << argv[2] << "\n";
    return -1;
  }
  std::cout << "Opening output file = " << output_file_name << "\n";

  float ratio = ((float)input_rate) / ((float)output_rate);
  int output_frame_samples =
      (int)ceil(ratio * num_channels * (float)input_frame_samples);
  std::vector<int16_t> input_data(input_frame_samples, 0);
  std::vector<int16_t> output_data(output_frame_samples, 0);
  int resampler_status = RESAMPLER_ERR_SUCCESS;

  // Create Resampler
  print_resampler_params(num_channels, input_rate, output_rate,
                         resampler_quality);
  std::unique_ptr<SpeexResamplerState, decltype(&speex_resampler_destroy)>
      resampler(speex_resampler_init(num_channels, input_rate, output_rate,
                                     resampler_quality, &resampler_status),
                speex_resampler_destroy);

  if (resampler_status != RESAMPLER_ERR_SUCCESS) {
    std::cerr << "Could not create SPEEX Resampler. Err = " << resampler_status;
    return -1;
  }

  // Read data from input_file, resample, write to output_file
  std::streamsize total_bytes_read = 0;
  std::cout << "Beginning to resample\n";
  while (true) {
    uint32_t input_samples = 0, num_processed_samples = 0;
    input_file.read(reinterpret_cast<char *>(input_data.data()),
                    input_frame_samples * sizeof(int16_t));

    if (input_file.gcount() % sizeof(int16_t)) {
      std::cerr << "Error reading from PCM file. Data is not aligned with "
                   "int16_t type\n";
      break;
    }
    if (input_file.gcount() == 0) {
      std::cout << "Reached EOF. Exiting event loop\n";
      break;
    }
    input_samples =
        static_cast<uint32_t>(input_file.gcount() / sizeof(int16_t));

    num_processed_samples = output_data.size();
    resampler_status = speex_resampler_process_interleaved_int(
        resampler.get(), input_data.data(), &input_samples, output_data.data(),
        &num_processed_samples);
    if (resampler_status != RESAMPLER_ERR_SUCCESS) {
      std::cerr << "Resampler Error: " << resampler_status << "\n";
      std::cerr << "Exiting\n";
      break;
    }

    output_file.write(reinterpret_cast<char *>(output_data.data()),
                      num_processed_samples * sizeof(int16_t));
    if (!input_file || !output_file) {
      std::cout << "Input and/or output file are not valid. Exiting\n";
      break;
    }
  }

  if (input_file) input_file.close();
  if (output_file) output_file.close();
  return 0;
}
