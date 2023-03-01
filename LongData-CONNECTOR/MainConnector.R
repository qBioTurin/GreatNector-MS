# library(devtools)
# install_github("qBioTurin/connector", ref="master",dependencies=TRUE)
library(connector)

rmarkdown::render("report.Rmd",
									output_dir = "Report",
									output_file ="EDSS.html",
									params = list(Data="EDSS",
																p = 4,
																g = c(3,4) ) )

rmarkdown::render("report.Rmd",
									output_dir = "Report",
									output_file =paste0("Treg_prod",".html"),
									params = list(Data="Treg_prod",
																p = 5,
																g = c(3,4)) )

rmarkdown::render("report.Rmd",
									output_dir = "Report",
									output_file =paste0("Th1_prod",".html"),
									params = list(Data="Th1_prod",
																p = 4,
																g = c(3:4)) )
rmarkdown::render("report.Rmd",
									output_dir = "Report",
									output_file =paste0("Th17_prod",".html"),
									params = list(Data="Th17_prod",
																p = 5,
																g = c(3:5)) )


rmarkdown::render("report.Rmd",
									output_dir = "Report",
									output_file =paste0("CD4",".html"),
									params = list(Data="CD4",
																	p = 5,
																	g = 3:4) )
