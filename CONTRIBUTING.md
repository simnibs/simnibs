Basics
======

Development takes place in the **develop** branch. Right after cloning, please do
```
git checkout develop
```

For a smoother experience, I've set-up a few hooks to help development. Please ensure you are using those by

```
git config core.hooksPath .githooks
```


Doing a New Release
====================


1. Bring changes from the Master branch by merging
```
git merge master
```
2. Checkout the master branch
``` 
git checkout master
```
3. Bring over changes from the develop branch using **checkout** **do NOT merge**
```
git checkout develop -- .
```
 This will stage everything, including the internal files. You can remove those with
```
git rm -r -f **/internal
```
During this process, some errors might occur, just remove the wrong files manually and proceed
4. Push to public and to origin
```
git push origin
git push -u public
```
