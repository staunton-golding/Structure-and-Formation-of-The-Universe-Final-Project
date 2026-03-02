# Structure and Formation of The Universe Final Project
Vanderbilt ASTR3800 final project (Spring 2021, Sophomore year)

This project entailed writing a friends-of-friends clustering algorithm for finding halos (groups of galaxies connected by a linking radius). 
This style of clustering algorithm differs from KNN / K-means, as it allows for trails of points with points in a halo potentially far closer to another
halo's centroid than to its own halo's centroid. The purpose of this algorithm is to reveal large-scale cosmic structures. It works by creating chains 
of points, all within a linking radius to one point already in the halo, as opposed to point clouds (where are points would be within a certain radius to its centroid.)
This project had some fun quirks - notably, the dataset was a "pac-man" style dataset, where the galaxies were in a toroidal space -- axes were assumed to be continuous (that is
the last point along the x-axis would be adjacent the first point on the x-axis). 

I am publishing this simply as a work example. Unfortunately, as this (school) project is from years ago, I could not source the .dat data file used. Apologies.

I included both the original livescript (.mlx), and a normal script (.m). MATLAB (online or application) is needed to view the .mlx file.
