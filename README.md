# cp2k developers at mcgill. branched from the official svn repository at r17181 (with minor ALMO EDA updates)

Purpose of creating this git repository is to keep track of numerous CP2K development projects at McGill.

Typical development cycle:

1. Master branch of this git repository is periodically syncronized with the trunk of the official svn repository.
2. Develop your own git branch, do not keep long-term incomplete projects in master git branch.
3. All git branches receive updates from (synchronized) master git branch.
4. When the development in a git branch is finished this branch is merged into master git branch.
5. All new developments in master git branch are send to the official svn repository.

================ Sychronization ===============
* goto dir that contains both .svn and .git

cd CP2K_HOME

* make sure you are on master branch

git checkout master

* make sure you master branch does not contain any uncomitted developments as indicated above it important not to keep long-term developments in master branch

git status

* compare your local svn to the official repository

svn status --show-updates

* get updates from the official repository

svn update

* check what changes git sees

git status

* decide which svn updates git master branch is allowed to receive, ideally all changes (except for obvious garbage). to stage all modified/deleted files use:

git add -u

* commit staged files

git commit -m "message here"

* send to github

git push origin master

==============================================
