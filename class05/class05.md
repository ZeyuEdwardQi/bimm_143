# Class 5: Data viz with ggplot2
Zeyu Qi(A17342618)

- [Background](#background)
- [Add some custom features](#add-some-custom-features)
- [Gene Expression](#gene-expression)
- [Going Further](#going-further)

## Background

There are many graphics systems in R for making plots and figures. These
include so-called *“base R” graphics* like the `plot()` function and add
on packages like **ggplot2**.

Let’s compare how we make a simple figure with these two systems:

We can use the in-built `cars` dataset:

``` r
head(cars)
```

      speed dist
    1     4    2
    2     4   10
    3     7    4
    4     7   22
    5     8   16
    6     9   10

``` r
plot(cars)
```

![](class05_files/figure-commonmark/unnamed-chunk-2-1.png)

Before I can use ggplot2 I need to install it on my computer. To do this
we can use the function `install.packages("ggplot2")`

> **N.B.** We never run `install.packages()` in our run quarto doc (we
> run it once only in our R console) as it would re-install the package
> every time we render.

Once installed we need to load up the package into our R brain:

``` r
library(ggplot2)
```

The main function in the **ggplot2** packages is called `ggplot()`

``` r
ggplot(cars)
```

![](class05_files/figure-commonmark/unnamed-chunk-4-1.png)

Every ggplot has at least 3 layers:

- the **data** (a datafrome)
- the **aes**thetics (how data maps to plot)
- the **geom** layer (how you want the plot drawn)

``` r
ggplot(cars) + 
  aes(x=speed, y=dist) +
  geom_point()
```

![](class05_files/figure-commonmark/unnamed-chunk-5-1.png)

## Add some custom features

Let’s add a trend line taht shows the relationship between speed and
distance.

``` r
ggplot(cars) + 
  aes(x=speed, y=dist) +
  geom_point() + 
  geom_smooth() +
  theme_bw() +
  labs(title="Stoping Dist of old cars", 
       x="Speed (MPH)",
       y="Distance (ft)")
```

    `geom_smooth()` using method = 'loess' and formula = 'y ~ x'

![](class05_files/figure-commonmark/unnamed-chunk-6-1.png)

> Q. Can you make the `geom_smoothe()` function produce a linear
> straight line fit to the data and turn off grey region.

------------------------------------------------------------------------

## Gene Expression

``` r
url <- "https://bioboot.github.io/bimm143_S20/class-material/up_down_expression.txt"
genes <- read.delim(url)
head(genes)
```

            Gene Condition1 Condition2      State
    1      A4GNT -3.6808610 -3.4401355 unchanging
    2       AAAS  4.5479580  4.3864126 unchanging
    3      AASDH  3.7190695  3.4787276 unchanging
    4       AATF  5.0784720  5.0151916 unchanging
    5       AATK  0.4711421  0.5598642 unchanging
    6 AB015752.4 -3.6808610 -3.5921390 unchanging

``` r
sum(genes$State == "up")
```

    [1] 127

A useful new function in this context is the `table()` function:

``` r
table(genes$State)
```


          down unchanging         up 
            72       4997        127 

``` r
ggplot(genes) +
  aes(Condition1, Condition2, col=State) +
  geom_point() + 
  scale_color_manual(values=c("red", "gray", "blue")) +
  theme_bw() +
  labs(x="No Drug", y="Drug", title="Expresion changes upon GLP-1 inhibitor treatment")
```

![](class05_files/figure-commonmark/unnamed-chunk-10-1.png)

## Going Further

Here we read the famous gapmider dataset:

``` r
# File location online
url <- "https://raw.githubusercontent.com/jennybc/gapminder/master/inst/extdata/gapminder.tsv"

gapminder <- read.delim(url)
head(gapminder)
```

          country continent year lifeExp      pop gdpPercap
    1 Afghanistan      Asia 1952  28.801  8425333  779.4453
    2 Afghanistan      Asia 1957  30.332  9240934  820.8530
    3 Afghanistan      Asia 1962  31.997 10267083  853.1007
    4 Afghanistan      Asia 1967  34.020 11537966  836.1971
    5 Afghanistan      Asia 1972  36.088 13079460  739.9811
    6 Afghanistan      Asia 1977  38.438 14880372  786.1134

> Q How many entries (i.e rows) are in this dataset?

``` r
nrow(gapminder)
```

    [1] 1704

> Q. How many different “country” entries are in this dataset?

``` r
length(table(gapminder$country))
```

    [1] 142

``` r
length(unique(gapminder$country))
```

    [1] 142

Let’s make our first plot of the entire dataset:

``` r
p <- ggplot(gapminder) + 
  aes(gdpPercap, lifeExp, col=continent) +
  geom_point(alpha=0.3) 
```

``` r
p
```

![](class05_files/figure-commonmark/unnamed-chunk-16-1.png)

I can add more layers to `p`

``` r
p + 
  facet_wrap(~continent)
```

![](class05_files/figure-commonmark/unnamed-chunk-17-1.png)

Make a plot for 1977 and 2007 only (not all the years in the dataset).

> Q. First use the **dplyr** package and the `filter()` function from
> that pacakge to extract the year 2007.

``` r
library(dplyr)
```

``` r
ggplot(filter(gapminder, year == 2007 | year == 1977)) + 
  aes(gdpPercap, lifeExp, col=continent, size=pop) +
  geom_point() + 
  facet_wrap(~year)
```

![](class05_files/figure-commonmark/unnamed-chunk-19-1.png)

> Q. Make a histogram of lifeExp faceted by continent

``` r
ggplot(gapminder) +
  aes(lifeExp) +
  geom_histogram() + 
  facet_wrap(~continent)
```

    `stat_bin()` using `bins = 30`. Pick better value `binwidth`.

![](class05_files/figure-commonmark/unnamed-chunk-20-1.png)

> Q. Make a histogram of lifeExp colored by continent

``` r
ggplot(gapminder) +
  aes(lifeExp, fill=continent) +
  geom_histogram()
```

    `stat_bin()` using `bins = 30`. Pick better value `binwidth`.

![](class05_files/figure-commonmark/unnamed-chunk-21-1.png)
