---
title: "Renaming camera trap images using Windows batch files"
date: 2022-06-06
categories:
  - blog
tags:
  - camera traps
  - tutorial
---

In this post, I explain how to add a prefix to and remove a suffix from multiple camera trap images using Windows batch files.

## Introduction

When working with hundreds of thousands of camera trap images, even simple things like renaming images can suddenly become a daunting task!

I encounter this situation pretty often during [my research](https://www.zooniverse.org/projects/peter-dot-stewart/prickly-pear-project-kenya), so I make use of Windows batch files to rename large numbers of images at once. In this post, I share how to do this.

Please note that **this will only work if you're on Windows**!

## VERY IMPORTANT NOTE

**Always work on a copy of your data. Never alter the original images, and keep the originals in a safe place away from where you are working, so that you have a backup in case things go wrong.**

**The code in this post worked for me, but I cannot guarantee that it will work for you - use it at your own risk. I always recommend testing batch files on a small subset of your images first, to confirm that the result is what you expected.**

## Adding a common prefix

When camera trap images come off the SD card, they usually have generic file names like `IMG_0001.jpg`. I often want to rename the images to something more useful by including the camera trap's site ID in the filename, for instance `Site_18_part2_IMG_0001.jpg`

To do this, I use this code:

```batch
FOR /f "delims=" %%F IN ('DIR /a-d /b *.jpg')  DO (
RENAME "%%F" "Site_18_part2_%%F")
```

This will append all of our images with the prefix `Site_18_part2_`. To change this prefix just alter the second line of code, making sure to keep the `%%F` bit.

To make this run, we have to create a `.bat` file in the folder which contains the camera trap images. **This file is going to rename every .jpg in the folder by adding the site prefix, so make sure that you move any images you don't want renamed to a separate folder!**

The first step is to create a `.txt` file in our folder. To do this, right click > New > Text Document. Name it something sensible - apparently it is important to avoid calling your batch file the same thing as other common batch files which may exist on your computer. In this example I've called my file`Prefix_camera_images.bat`

We then paste the code above into this `.txt` file, and save it.

We then need to turn our `.txt` file into a batch file. We can do this by using save as - we need to select "All files" in the drop-down menu, and type `.bat` at the end of the file name:

![](assets/images/post_images/useful_batch_files/save_as.jpg)

Then we save the file. You'll see that the batch file has now appeared in the folder with our camera trap images:

![](assets/images/post_images/useful_batch_files/prefix_before.jpg)

Double-click it, and all of your images will be renamed to include the prefix:

![](assets/images/post_images/useful_batch_files/prefix_after.jpg)

If you've got a large number of images this might take a few seconds - but it's still MUCH faster than doing it manually!

After it's done, it's a good idea to **delete the batch file** as if you accidentally double-click it again, it's going to add an unwanted second prefix to your images!

## Removing a bracketed number

I recently encountered a situation in which images from different camera trap sites had been placed in the same folder - this resulted in many of the images being appended with a space followed by a bracketed number, for example:

![](assets/images/post_images/useful_batch_files/before.jpg)

I wanted to get rid of the space and bracketed number, so that e.g. `IMG_0001 (3).jpg` would become `IMG_0001.jpg`

To do this, I used the following code:

```batch
setlocal enableDelayedExpansion
FOR /f "delims=" %%F IN ('DIR /a-d /b *.jpg') DO (
  set filename="%%~nxF"
  ren "%%F" "!filename: (3)=!"
)
```

If you want to remove a different suffix, change the bit on the fourth line between `!filename:` and `=!` - note that there is a space between `!filename:` and `(3)` in my example, because I wanted to remove the space before the bracketed number as well.

We run this in the same way as for the prefix batch file above - check that section of the post for the step-by-step instructions for creating and running the batch file.

Here is the result:

![](assets/images/post_images/useful_batch_files/after.jpg)
