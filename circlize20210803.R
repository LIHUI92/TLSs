library(circlize)

dat <- read.csv('C:/Users/lihui/Desktop/GEEA_result.csv', header = T, sep = ',',row.names = 1)

rownames(dat) <- dat$ID

#以category排序

dat <- dat[order(dat$category),]
circos.par(gap.degree = 0.5, start.degree = 90)
plot_data <- dat[c('ID', 'gene_num.min', 'gene_num.max')]

#分配颜色

color_assign <- c('A. Immune receptors' = '#954572', 'Cytokines' = '#F7CC13', 'B. Others'='#33A02CFF')
ko_color <- color_assign[dat$category]

circos.genomicInitialize(plot_data, plotType = NULL, major.by = 1)
circos.track(
  ylim = c(0, 1), track.height = 0.05, bg.border = NA, bg.col = ko_color,
  panel.fun = function(x, y) {
    ylim = get.cell.meta.data('ycenter')#ylim、xlim
    xlim = get.cell.meta.data('xcenter')
    sector.name = get.cell.meta.data('sector.index')#sector.name
    circos.axis(h = 'top', labels.cex = 0.4, labels.niceFacing = FALSE)
    circos.text(xlim, ylim, sector.name, cex = 0.4, niceFacing = FALSE)
  } )

##第二圈，绘制富集的基因和富集p值
plot_data <- dat[c('ID', 'gene_num.min', 'Count', 'log10Pvalue')]
label_data <- dat['Count']
p_max <- round(max(dat$'log10Pvalue')) + 1  
colorsChoice <- colorRampPalette(c('#FF906F', '#861D30'))
color_assign <- colorRamp2(breaks = 0:p_max, col = colorsChoice(p_max + 1))

circos.genomicTrackPlotRegion(
  plot_data, track.height = 0.08, bg.border = NA, stack = TRUE,
  panel.fun = function(region, value,...) {
    circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA,...)
    ylim = get.cell.meta.data('ycenter')  
    xlim = label_data[get.cell.meta.data('sector.index'),1] / 2
    sector.name = label_data[get.cell.meta.data('sector.index'),1]
    #circos.text(xlim, ylim, sector.name, cex = 0.4, niceFacing = FALSE)
  } )

##第三圈，绘制上下调基因
dat$all.regulated <- dat$up.regulated + dat$down.regulated
dat$up.proportion <- dat$up.regulated / dat$all.regulated
dat$down.proportion <- dat$down.regulated / dat$all.regulated

dat$up <- dat$up.proportion * dat$gene_num.max
plot_data_up <- dat[c('ID', 'gene_num.min', 'up')]
names(plot_data_up) <- c('ID', 'start', 'end')
plot_data_up$type <- 1 

dat$down <- dat$down.proportion * dat$gene_num.max + dat$up
plot_data_down <- dat[c('ID', 'up', 'down')]
names(plot_data_down) <- c('ID', 'start', 'end')
plot_data_down$type <- 2 

plot_data <- rbind(plot_data_up, plot_data_down)
label_data <- dat[c('up', 'down', 'up.regulated', 'down.regulated')]
color_assign <- colorRamp2(breaks = c(1, 2), col = c('red', 'blue'))

circos.genomicTrackPlotRegion(
  plot_data, track.height = 0.08, bg.border = NA, stack = TRUE, 
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...) 
    ylim = get.cell.meta.data('cell.bottom.radius') - 0.5 
    xlim = label_data[get.cell.meta.data('sector.index'),1] / 2
    sector.name = label_data[get.cell.meta.data('sector.index'),3]
    #circos.text(xlim, ylim, sector.name, cex = 0.4, niceFacing = FALSE)  
    xlim = (label_data[get.cell.meta.data('sector.index'),2]+label_data[get.cell.meta.data('sector.index'),1]) / 2
    sector.name = label_data[get.cell.meta.data('sector.index'),4]
    #circos.text(xlim, ylim, sector.name, cex = 0.4, niceFacing = FALSE)
  } )

##第四圈，绘制富集因子
plot_data <- dat[c('ID', 'gene_num.min', 'gene_num.max', 'rich.factor')] 
label_data <- dat['category']  
color_assign <- c('A. Immune receptors' = '#954572', 'Cytokines' = '#F7CC13', 'B. Others'='#33A02CFF')#各二级分类的名称和颜色

#value值超出了track.height，乘以该值，控制bar的高度即可

circos.genomicTrack(
  plot_data, ylim = c(0, 1), track.height = 0.3, bg.col = 'gray95', bg.border = NA,
  panel.fun = function(region, value, ...) {
    sector.name = get.cell.meta.data('sector.index')  #sector.name 
    circos.genomicRect(region, value*0.3, col = color_assign[label_data[sector.name,1]], border = NA, ytop.column = 1, ybottom = 0, ...) 
    circos.lines(c(0, max(region)), c(0.5, 0.5), col = 'gray', lwd = 0.3) 
  } )

#画legend需要这个包（注意替换自己的变量）

library(ComplexHeatmap)

category_legend <- Legend(
  labels = c('A. Immune receptors', 'Cytokines', 'B. Others'),#各二级分类的名称
  type = 'points', pch = NA, background = c('#954572', '#F7CC13', '#33A02CFF'), #各二级分类的颜色
  labels_gp = gpar(fontsize = 8), grid_height = unit(0.5, 'cm'), grid_width = unit(0.5, 'cm'))

updown_legend <- Legend(
  labels = c('Up-regulated', 'Down-regulated'), 
  type = 'points', pch = NA, background = c('red', 'blue'), 
  labels_gp = gpar(fontsize = 8), grid_height = unit(0.5, 'cm'), grid_width = unit(0.5, 'cm'))

pvalue_legend <- Legend(
  col_fun = colorRamp2(round(seq(0, p_max, length.out = 6), 0), 
                       colorRampPalette(c('#FF906F', '#861D30'))(6)),
  legend_height = unit(3, 'cm'), labels_gp = gpar(fontsize = 8), 
  title_gp = gpar(fontsize = 9), title_position = 'topleft', title = '-Log10(Pvalue)')

lgd_list_vertical <- packLegend(category_legend, updown_legend, pvalue_legend)
grid.draw(lgd_list_vertical)

circos.clear()
dev.off()