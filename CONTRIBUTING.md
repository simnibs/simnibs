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


1. Make sure that the build is working

2. Checkout the release branch
``` 
git checkout release
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

4. Push to public
```
git push -u public
```
