# Pull Request Template for OpenMD

**IMPORTANT: Please do not create a Pull Request without creating an issue first.**

*Any change needs to be discussed before proceeding. Failure to do so may result in the rejection of the pull request.*

Please provide enough information so that others can review your pull request. You can skip this if you're just fixing a typo.

Explain the **details** for making this change. What existing problem does the pull request solve?

Ex:

1. If you "Added a/changed the function to do X", explain why:

    - It is necessary to have a way to do X.
    - If there already exists a way, why is your implementation better.

2. If you "Fixed bug/error in X", explain:

    - What was the bug/error (if you already made an issue, please link to it here).
    - How does your implementation fix the issue.

## Code style and formatting

Attempt to match the surrounding code styling as much as possible. Check the [Contributors Style Guidelines section](https://github.com/OpenMD/OpenMD/blob/main/.github/CONTRIBUTING.md#Style-guidelines) for how to write your code and the [Contributors Code Formatting section](https://github.com/OpenMD/OpenMD/blob/main/.github/CONTRIBUTING.md#Code-formatting) for how to format your code.

### Closing Issues

Put `closes #XXXX` in your comment to auto-close the issue that your PR fixes (if such).

---

## Proposed changes

-
-
-

## Motivation behind changes

### Test plan

Demonstrate the code is solid. Example: The exact commands you ran and their output, screenshots / videos if the pull request changes UI.

*Make sure tests pass on all of the [relevant CI workflows](https://github.com/OpenMD/OpenMD/blob/main/.github/workflows/build.yml).*

### Pull Request Readiness Checklist

See details at [CONTRIBUTING.md](https://github.com/OpenMD/OpenMD/blob/main/.github/CONTRIBUTING.md).

- [ ] I agree to contribute to the project under OpenMD's [BSD 3-Clause License](https://github.com/OpenMD/OpenMD/blob/main/LICENSE).

- [ ] To the best of my knowledge, the proposed patch is not based on a code under GPL or other license that is incompatible with OpenMD.

- [ ] The PR is proposed to proper branch.

- [ ] There is reference to original bug report and related work.

- [ ] There is accuracy test, performance test and test data in the repository, if applicable.

- [ ] The feature is well documented and sample code can be built with the project CMake.
