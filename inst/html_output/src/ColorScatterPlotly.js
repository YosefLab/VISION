/*
   Initializes a zoomable scatter plot in the element "parent"
   parent = , for example, "#chart_div"
*/
function ColorScatter(parent, colorbar, legend)
{
    this.parent = $(parent).get(0)
    this.width = $(parent).width()
    this.height = $(parent).height()
}

ColorScatter.prototype.clearData = function() {
}

/*
 * points: array of [x, y, color, name]
 *
 */
ColorScatter.prototype.setData = function(points, isFactor,
    full_color_range, selected_cells, diverging_colormap)
{

    var circle_radius = this.pointsToRadius(points.length)

    var data = []; // Holds plotly traces
    if(isFactor){

        for (var j =0; j < 4; j+= 1){
            var X = [],
                Y = [],
                id_sub = [],
                n = 1000,
                i;

            for (i = 0; i < n; i += 1) {
                X.push(points[i+n*j][0])
                Y.push(points[i+n*j][1])
                id_sub.push('ok')
            }

            var trace = {
                x: X,
                y: Y,
                mode: "markers",
                type: "scattergl",
                text: id_sub,
                name: j.toString(),
                marker: {
                    size: 4.2343,
                }
            }

            data.push(trace)
        }

    } else if(isFactor) {
        var c = _.map(points, p => p[2])
        var unique = d3.set(c).values();
        _.forEach(unique, level => {
            var subset = _.filter(points, p => p[2] == level)
            var x_sub = _.map(subset, p => p[0])
            var y_sub = _.map(subset, p => p[1])
            var id_sub = _.map(subset, p => p[3])
            var trace = {
                x: x_sub,
                y: y_sub,
                mode: 'markers',
                type: 'scattergl',
                text: id_sub,
                name: level.toString(),
                marker: {
                    size: circle_radius,
                },
            }
            data.push(trace)
        });

    } else {

        var x = _.map(points, p => p[0])
        var y = _.map(points, p => p[1])
        var c = _.map(points, p => p[2])
        var id = _.map(points, p => p[3])

        var cvals = c.filter(cv => cv !== "NA")
        cvals.sort(d3.ascending); // Needed for quantile
        var low, high;
        if(full_color_range){
            low = d3.min(cvals)
            high = d3.max(cvals)
        } else {
            low = d3.quantile(cvals, 0.02)
            high = d3.quantile(cvals, 0.98)
        }

        var colorscale
        if(diverging_colormap){
            colorscale = 'Viridis'
        } else {
            colorscale = 'Greys'
        }

        var trace1 = {
            x: x,
            y: y,
            mode: 'markers',
            type: 'scattergl',
            text: id,
            marker: {
                size: circle_radius,
                color: c,
                cmin: low,
                cmax: high,
                colorbar: {
                    title: 'Colorbar'
                },
                colorscale: colorscale,
            },

        }

        data.push(trace1)

    }

    Plotly.react(this.parent, data)
}

ColorScatter.prototype.pointsToRadius = function(n_points)
{
    // Pick a point size based on the number of scatter
    // plot points

    // Using a piecewise function for this
    var a, b


    var scale_factor = Math.min(this.width, this.height)/600

    if(n_points <= 10000){
        a = 116500
        b = .284
    } else {
        a = 19800
        b = .748
    }

    return Math.pow(a/n_points, b) * scale_factor

}

ColorScatter.prototype.setTreeData = function(tree_points, tree_adj)
{

    this.tree_points = tree_points
    this.tree_adj = tree_adj
}

/*
 * Automatically adjusts the scale if the points are too wide
 * or too narrow
 */
ColorScatter.prototype.autoZoom = function() {

};



ColorScatter.prototype.hover_cells = function(cell_ids) {

};



ColorScatter.prototype.getSelected = function() {

}
