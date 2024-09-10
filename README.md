## Research Software Engineer

Developing Teal Apps for interactive data analysis.

`{r} webr::install("teal")`

```{r}
install.packages("shinylive")
shinylive::export("myapp", "site")
httpuv::runStaticServer("site")
```
