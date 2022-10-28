avg.nodes <- function (object)
{
num <- nrow(x = object@meta.data)
res0 <- object@meta.data$integrated_snn_res.0
res0.1 <- object@meta.data$integrated_snn_res.0.1
res0.2 <- object@meta.data$integrated_snn_res.0.2
res0.3 <- object@meta.data$integrated_snn_res.0.3
res0.4 <- object@meta.data$integrated_snn_res.0.4
res0.5 <- object@meta.data$integrated_snn_res.0.5
res0.6 <- object@meta.data$integrated_snn_res.0.6
res0.7 <- object@meta.data$integrated_snn_res.0.7
res0.8 <- object@meta.data$integrated_snn_res.0.8
res0.9 <- object@meta.data$integrated_snn_res.0.9
res1 <- object@meta.data$integrated_snn_res.1
res1.1 <- object@meta.data$integrated_snn_res.1.1
res1.2 <- object@meta.data$integrated_snn_res.1.2
res1.3 <- object@meta.data$integrated_snn_res.1.3
res1.4 <- object@meta.data$integrated_snn_res.1.4
res1.5 <- object@meta.data$integrated_snn_res.1.5

mat <- matrix(c(res0,res0.1,res0.2,res0.3,res0.4,res0.5,res0.6,res0.7,res0.8,
res0.9,res1,res1.1,res1.2,res1.3,res1.4,res1.5), nrow = num, dimnames = list(c(rownames(object@meta.data)),
c("res0","res0.1","res0.2","res0.3","res0.4","res0.5","res0.6","res0.7","res0.8","res0.9","res1",
"res1.1","res1.2","res1.3","res1.4","res1.5")))

nodes <- get_tree_nodes(clusterings = mat, prefix = "res", node_aes_list = FALSE, metadata = FALSE)
nodes.agg <- aggregate(nodes$sc3_stability,by=list(res=nodes$res),data=nodes,FUN=mean)
return(nodes.agg)
}
