# Other `BioFormatsImage` getters

`isFULL` indicates if the extent is the full extent of the image.
`origin` gets the x-y coordinates of the origin of the image, i.e. the
smallest possible x-y coordinate values within the full image.

## Usage

``` r
# S4 method for class 'BioFormatsImage'
isFull(x)

# S4 method for class 'BioFormatsImage'
origin(x)

# S4 method for class 'BioFormatsImage'
transformation(x)
```

## Arguments

- x:

  A
  [`BioFormatsImage`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/BioFormatsImage.md)
  object.

## Value

For `isFull`: Logical scalar indicating whether the extent is the full
extent. For `origin`: Numeric vector of length 2. For `transformation`,
a list.
