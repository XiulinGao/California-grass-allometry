library(LeafArea)



home_path = path.expand("~")
img_path  = file.path(home_path,"Downloads/allomery-leaf-area-scan/")
out_path  = file.path(home_path,"California-grass-allometry/data/")
lfarea    = run.ij(set.directory  = img_path
                  ,distance.pixel = 826
                  ,known.distance  = 21
                  ,log=TRUE)
mean_larea = lfarea$summary
each_lf = tibble::as_tibble(lfarea$each.image)
each_area = tibble::as_tibble(unlist(lfarea$each.image))
nam = names(each_lf)
nam = gsub(pattern=".jpeg.txt",replacement="",x=nam)
each_area$sample = nam
f_output = "leaf-area.csv"
fwrite(each_area,file = file.path(out_path,f_output))
