## Repository
Development takes place on our [Gitlab Repository](https://gitlab.xiph.org/xiph/speexdsp). You will need to register an account first and get approval from an admin to be able to fork the repository. If you do not get approval in a timely manner ping on IRC at #xiph or #speex @ irc.libera.chat.

There is also a [Github mirror](https://github.xiph.org/xiph/speexdsp) for visibility/convenience but **pull requests are discouraged**, as it means maintainers will have to recreate them as merge requests on our gitlab to ensure they pass CI.

## Development

Once you have forked the repository, you can submit merge requests (using the branch from your repository as the source branch and the `main` branch of our repository as the target branch).

We expect all commits to be GPG-signed (see https://docs.gitlab.com/user/project/repository/signed_commits/gpg/). We recommend making this setting the default for this repo using `git config --local commit.gpgsign true`.

Once a merge request passes CI (including verifying that all the commits are signed) and reviewed, it can be merged. If no one has responded in 48 hours, ping on IRC at #xiph or #speex @ irc.libera.chat

We also accept patches if Gitlab proves too cumbersome (although again, merge requests are preferred).

## Bug reporting

Issues can be filed using the [Work Items](https://gitlab.xiph.org/xiph/speexdsp/-/work_items) dialogue by choosing **New Item** and **Issue**.

**Security bugs should be filed as confidential issues.** There's a checkbox to this effect at the bottom of the **Issue** form. 

## Project History

Originally **SpeexDSP** (the DSP library) and **Speex** (the codec library) were a single project. They've since been split. For **Speex** codec related development, see [its Gitlab](https://gitlab.xiph.org/speex).
