make_layout <- function(m){
    start_layouts <- list()
    start_layout_inters <- 1
    n <- 1

    while (start_layout_inters[length(start_layout_inters)] != 0 & n < 10) {
        net <- intergraph::asNetwork(m)
        xy <- network::as.matrix.network.adjacency(net)
        layout1 <- sna::gplot.layout.kamadakawai(xy, layout.par=list(niter=500))
        # layout1 <- layout_with_kk(m) # if (layout=="igraph_kk")

        start_layout_inters <- c(start_layout_inters,
                                 n_intersect_segm(store_all_info(m, layout1)$lines_to_check))

        start_layouts[[n]] <- layout1
        n <- n + 1
    }

    layout1 <- start_layouts[[which.min(start_layout_inters[-1])]]
}

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

store_all_info <- function(m, layout1){

    nt <- as_data_frame(m, "vertices")
    nt$x <- layout1[,1]
    nt$y <- layout1[,2]

    et <- as_data_frame(m, "edges")

    et$x <- sapply(seq(et$from), function(i) nt$x[nt$name==et$from[i]])
    et$y <- sapply(seq(et$from), function(i) nt$y[nt$name==et$from[i]])
    et$s_num <- sapply(seq(et$from), function(i) which(nt$name==et$from[i]))

    et$xend <- sapply(seq(et$to), function(i) nt$x[nt$name==et$to[i]])
    et$yend <- sapply(seq(et$to), function(i) nt$y[nt$name==et$to[i]])
    et$e_num <- sapply(seq(et$to), function(i) which(nt$name==et$to[i]))

    lines_to_check <- data.frame(
        x=et$x,
        y=et$y,
        xend=et$xend,
        yend=et$yend,
        IDstart=as.character(et$from),
        IDend=as.character(et$to),
        IDs_num=et$s_num,
        IDe_num=et$e_num)

    lines_to_check$middle_x <- (lines_to_check$x + lines_to_check$xend) / 2
    lines_to_check$middle_y <- (lines_to_check$y + lines_to_check$yend) / 2

    return(list(
        nt=nt,
        et=et,
        lines_to_check=lines_to_check)
    )
}

# segment intersection
on_segment <- function(x, y, xend, yend, rx, ry){
    return (xend <= max(x, rx) & xend >= min(x, rx) & yend <= max(y, ry) & yend >= min(y, ry))
}

orientation <- function(x, y, xend, yend, rx, ry){
    sign((yend - y) * (rx - xend) - (xend - x) * (ry - yend))
}

intersect_segm <- function(segm1, segm2){
    o1 <- orientation(segm1$x, segm1$y, segm1$xend, segm1$yend, segm2$x, segm2$y)
    o2 <- orientation(segm1$x, segm1$y, segm1$xend, segm1$yend, segm2$xend, segm2$yend)
    o3 <- orientation(segm2$x, segm2$y, segm2$xend, segm2$yend, segm1$x, segm1$y)
    o4 <- orientation(segm2$x, segm2$y, segm2$xend, segm2$yend, segm1$xend, segm1$yend)

    return(
        (o1!=o2 & o3!=o4) |
            (o1==0 & on_segment(segm1$x, segm1$y, segm2$x, segm2$y, segm1$xend, segm1$yend)) |
            (o2==0 & on_segment(segm1$x, segm1$y, segm2$xend, segm2$yend, segm1$xend, segm1$yend)) |
            (o3==0 & on_segment(segm2$x, segm2$y, segm1$x, segm1$y, segm2$xend, segm2$yend)) |
            (o4==0 & on_segment(segm2$x, segm2$y, segm1$xend, segm1$yend, segm2$xend, segm2$yend)))
}

n_intersect_segm <- function(lines_to_check){
    inter_matrix <- outer(seq_len(nrow(lines_to_check)), seq_len(nrow(lines_to_check)),
                          function(i, j){
                              ifelse((as.character(lines_to_check$IDstart[i])==as.character(lines_to_check$IDstart[j])) |
                                         (as.character(lines_to_check$IDend[i])==as.character(lines_to_check$IDend[j])) |
                                         (as.character(lines_to_check$IDstart[i])==as.character(lines_to_check$IDend[j])) |
                                         (as.character(lines_to_check$IDend[i])==as.character(lines_to_check$IDstart[j])),
                                     F, intersect_segm(lines_to_check[i, 1:4], lines_to_check[j, 1:4]))
                          })
    inter_matrix_ltri_vector <- inter_matrix[lower.tri(inter_matrix)]
    return(length(inter_matrix[inter_matrix_ltri_vector==T]))
}

# grobs params production
get_nlabel_boxes <- function(layout1, nlabel_semisizes){
    return(data.frame("x1"=layout1[, 1] - nlabel_semisizes$semi_w,
                      "y1"=layout1[, 2] - nlabel_semisizes$semi_h,
                      "x2"=layout1[, 1] + nlabel_semisizes$semi_w,
                      "y2"=layout1[, 2] + nlabel_semisizes$semi_h
    )
    )
}

get_elabel_boxes <- function(layout1, edges, elabel_semisizes){
    return(data.frame("x1"=(layout1[edges[,1], 1] + layout1[edges[,2], 1]) / 2 -
                          elabel_semisizes$semi_w,
                      "y1"=(layout1[edges[,1], 2] + layout1[edges[,2], 2]) / 2 -
                          elabel_semisizes$semi_h,
                      "x2"=(layout1[edges[,1], 1] + layout1[edges[,2], 1]) / 2 +
                          elabel_semisizes$semi_w,
                      "y2"=(layout1[edges[,1], 2] + layout1[edges[,2], 2]) / 2 +
                          elabel_semisizes$semi_h
    )
    )
}

# force algorithm utils
intersect_box <- function(box1, box2){
    return (box1$x1 < box2$x2 & box1$y1 < box2$y2 & box1$x2 > box2$x1 & box1$y2 > box2$y1)
}

quad_dist <- function(c1, c2){
    return ((c2$x - c1$x)^2 + (c2$y - c1$y)^2)
}

semiwidth <- function(box){
    return ((box$x2 - box$x1) / 2)
}

semiheight <- function(box){
    return((box$y2 - box$y1) / 2)
}

center <- function(c1, c2){
    return (data.frame("x"=(c1[1] + c2[1]) / 2, "y"=(c1[2] + c2[2]) / 2))
}

centroid <- function(box){
    return (data.frame("x"=box$x1 + semiwidth(box), "y"=box$y1 + semiheight(box)))
}

unit_vector <- function(c1, c2){
    return (data.frame("x"=(c1$x - c2$x) / sqrt(quad_dist(c1, c2)), "y"=(c1$y - c2$y) / sqrt(quad_dist(c1, c2))))
}


.messagef <- function (...)  { message(sprintf(...)) }
.warningf <- function (...)  { warning(sprintf(...)) }

.replaceNA <- function(x, y) { ifelse(is.na(x), y, x) }

.pairwiseCompare <- function (FUN, list1, list2 = list1, ...)
{
    additionalArguments <- list(...)
    f1 <- function(...) {
        mapply(FUN = function(x, y) {
            do.call(FUN, c(list(list1[[x]], list2[[y]]), additionalArguments))
        }, ...)
    }
    z <- outer(seq_along(list1), seq_along(list2), FUN = f1)
    rownames(z) <- names(list1)
    colnames(z) <- names(list2)
    z
}

.intersectionSize <- function(...) { length(intersect(...))}
