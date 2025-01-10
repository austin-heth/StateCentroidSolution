# StateCentroidSolution
Written in Python. This code implements my algorithm for obtaining the latitude and longitude coordinate corresponding to the centroid of a geographical surface -- in this case, the state of South Carolina.

# Design
The approach I took to solving this problem is summarized as follows:
My calculations rely on data obtained from the US Federal Communications Commission (https://www.fcc.gov/media/radio/us-county-overlays-kml) which provides a list of connecting coordinates (Lat, Long) associated with the boundaries of every county of a given state.
Using this data, I then determine the centroid and approx. surface area of each county. I do this by "pixelating" the county, i.e., deconstructing the 2-D space into little pixels of a pre-defined resolution (in this case, I choose a pixel size of 0.00001° lat. x 0.00001° long.).
Once the space is pixelated, I then identify which pixels fall within the boundary of the county and use the weighted average (adjusting the sq. km. measure of each pixel based on its level of latitude).
After I determine all the counties' centroids and surface areas, I compute the weighted average of these to obtain the centroid and approx. total surface area of the state.
