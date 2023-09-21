force_alg <- function(layout1,
                      nlabel_semisizes,
                      elabel_semisizes,
                      edges,
                      xlim, ylim,
                      n_iter, force){

    longer_goes_Y <- mean(nlabel_semisizes[ ,1])
    intersection <- TRUE
    iter <- 1

    layouts <- list()

    while(iter <= n_iter & intersection){

        # Let's try to add some random, firstly, to the whole layout1
        # (to get out of the local minima) # so called "jitter"
        # layout1[,1] <- layout1[,1] + 0.5 * runif(1, min=-mean(nlabel_semisizes$semi_w), max=mean(nlabel_semisizes$semi_w))
        # layout1[,2] <- layout1[,2] + 0.1 * runif(1, min=-mean(nlabel_semisizes$semi_w), max=mean(nlabel_semisizes$semi_w))

        intersection <- FALSE
        nlabel_boxes <- get_nlabel_boxes(layout1, nlabel_semisizes)
        elabel_boxes <- get_elabel_boxes(layout1, edges, elabel_semisizes)
        force_layout <- matrix(0, nrow=nrow(nlabel_boxes), ncol=2)

        # FORCE ACCUMULATION
        # NODES
        # ...and nodes
        inter_matrix1 <- outer(seq_len(nrow(nlabel_boxes)), seq_len(nrow(nlabel_boxes)),
                               function(i,j) { intersect_box(nlabel_boxes[i,], nlabel_boxes[j,]) })
        diag(inter_matrix1) <- FALSE
        inters <- which(inter_matrix1, arr.ind=TRUE)

        if(length(inters)!=0){
            intersection <- TRUE

            c1.1 <- centroid(nlabel_boxes[inters[,1],])
            c2.1 <- centroid(nlabel_boxes[inters[,2],])

            f.1 <- force * unit_vector(c1.1, c2.1) / pmax(quad_dist(c1.1, c2.1), 0.01)
            f.1 <- cbind(who=inters[,1], f.1)
            f.sum.1 <- aggregate(. ~ who, data = f.1, sum)
            force_layout[f.sum.1$who, ] <- force_layout[f.sum.1$who, ] + as.matrix(f.sum.1[,2:3])

            force_layout[, 2][nlabel_semisizes[,1] > longer_goes_Y] <- 2 * force_layout[, 2][nlabel_semisizes[,1] > longer_goes_Y]
        }

        # ...and edges
        inter_matrix2 <- outer(seq_len(nrow(nlabel_boxes)), seq_len(nrow(elabel_boxes)),
                               function(i,j) { intersect_box(nlabel_boxes[i,], elabel_boxes[j,]) })
        inters2 <- which(inter_matrix2, arr.ind=TRUE)

        if(length(inters2)!=0){
            intersection <- TRUE

            c1.2 <- centroid(nlabel_boxes[inters2[,1],])
            c2.2 <- centroid(elabel_boxes[inters2[,2],])

            f.2 <- force * unit_vector(c1.2, c2.2) / pmax(quad_dist(c1.2, c2.2), 0.01)
            f.2 <- cbind(who=inters2[,1], f.2)
            f.sum.2 <- aggregate(. ~ who, data=f.2, sum)
            force_layout[f.sum.2$who, ] <- force_layout[f.sum.2$who, ] + as.matrix(f.sum.2[,2:3])

            force_layout[, 2][nlabel_semisizes[,1] > longer_goes_Y] <- 2 * force_layout[, 2][nlabel_semisizes[,1] > longer_goes_Y]
        }

        # EDGES
        # ...and edges
        inter_matrix3 <- outer(seq_len(nrow(elabel_boxes)), seq_len(nrow(elabel_boxes)),
                               function(i,j) { intersect_box(elabel_boxes[i,], elabel_boxes[j,]) })
        diag(inter_matrix3) <- FALSE
        inters3 <- which(inter_matrix3, arr.ind=TRUE)

        if(length(inters3)!=0){
            intersection <- TRUE

            c1.3 <- centroid(elabel_boxes[inters3[,1], ])
            c2.3 <- centroid(elabel_boxes[inters3[,2], ])

            f.3 <- force * unit_vector(c1.3, c2.3) / pmax(quad_dist(c1.3, c2.3), 0.01)
            f.3 <- rbind(cbind(who=edges[inters3[,1], 1], f.3), cbind(who=edges[inters3[,1], 2], f.3))
            f.sum.3 <- aggregate(. ~ who, data = f.3, sum)
            force_layout[f.sum.3$who, ] <- force_layout[f.sum.3$who, ] + as.matrix(f.sum.3[,2:3])
        }

        # Prohibit strong displacements
        force_layout[force_layout > 0.1] <- 0.1
        force_layout[force_layout < -0.1] <- -0.1


        # FORCE REALIZATION
        if (iter!=n_iter){

            # Friction simulation
            layout1 <- layout1 + force_layout * (1 - (1e-3) * iter)

            # Prohibit crawling away borders
            layout1[,1][layout1[,1] < xlim] <- xlim
            layout1[,1][layout1[,1] > 1 - xlim] <- 1 - xlim
            layout1[,2][layout1[,2] < ylim] <- ylim
            layout1[,2][layout1[,2] > 1 - ylim] <- 1 - ylim

            layouts[[iter]] <- layout1
        }
        iter <- iter + 1
    }

    # data from (n-1) iteration:
    return(
        list(layout=layout1,
             force_ends=layout1 + force_layout,
             nlabel_boxes=nlabel_boxes,
             elabel_boxes=elabel_boxes,
             layouts=layouts
        ))
}
