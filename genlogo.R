#install.packages("hexSticker")
library(hexSticker)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("emojifont")
library(emojifont)
# fa4
#r <- ggplot() + geom_fontawesome("fa-sitemap", color='red', size = 20) + theme_void()

# fa5
#install.packages("fontawesome")
library(fontawesome)
aaa=fa(name = "fas fa-sitemap", fill = "#b1116b")


print_fontawesome<- function (x, view = interactive(), ...) 
{
  dots <- list(...)
  html <- paste(x, collapse = "\n")
  c("<!DOCTYPE html>", "<html>", "<head>", "<meta charset=\"utf-8\">", 
    "</head>", "<body>", html, "</body>", "</html>") %>% 
    paste(collapse = "\n") %>% htmltools::HTML() %>% htmltools::html_print()
  return(html)
}

library(rsvg)
cat(print_fontawesome(aaa, view=FALSE),file="man/figures/fa-sitemap.svg")
SVGres <- rsvg_svg("man/figures/fa-sitemap.svg","man/figures/Cairo-fa-sitemap.svg")
SVGres <- rsvg_svg("man/figures/fa-sitemap.svg")
PDFres <- rsvg_pdf("man/figures/fa-sitemap.svg","man/figures/Cairo-fa-sitemap.pdf")
PDFres <- rsvg_pdf("man/figures/fa-sitemap.svg")

library(png)
library(grid)
library(cowplot)
library(magick)
theme_set(theme_cowplot())

r = ggdraw() +
  draw_image(image_fill(
    image_fill(
      image_fill(
        image_fill(
          image_fill(
            image_fill(
              image_fill(
                image_fill(
                  image_fill(
                    image_fill(
                      image_fill(image_append(c(
                        image_rotate(image_background(
                          image_read_pdf("man/figures/Cairo-fa-sitemap.pdf"), 'black'
                        ), 90),
                        image_rotate(image_background(
                          image_crop(
                            image_read_pdf("man/figures/Cairo-fa-sitemap.pdf"),
                            "2667x1450+0+683"
                          ), 'black'
                        ), 270)
                      )),
                      "orange", point = "+100+200", fuzz = 20),
                      "purple",
                      point = "+3000+2500",
                      fuzz = 20
                    ),
                    "lightgreen",
                    point = "+100+2500",
                    fuzz = 20
                  ),
                  "blue",
                  point = "+3000+200",
                  fuzz = 20
                ),
                "#EE3B3B",
                point = "+3000+1350",
                fuzz = 20
              ),
              "#00C5CD",
              point = "+100+1350",
              fuzz = 20
            ),
            "none",
            point = "+1+1",
            fuzz = 20 
          ),
          "none",
          point = "+3582+2666",
          fuzz = 20
        ),
        "none",
        point = "+3582+1",
        fuzz = 20
      ),
      "none",
      point = "+1+2666",
      fuzz = 20
    ),
    "none",
    point = "+1330+1",
    fuzz = 20
  ),
  scale = .55)

sticker(
  r,
  package = "Patterns",
  p_size = 8,
  s_x = .98,
  s_y = 0.7,
  s_width = 1.7,
  s_height = 1.3,
  p_x = 1,
  p_y = 1.3,
  url = "https://cran.r-project.org/package=Patterns",
  u_color = "white",
  u_size = 1.1,
  h_fill = "black",
  h_color = "grey",
  filename = "man/figures/logo.png"
)
