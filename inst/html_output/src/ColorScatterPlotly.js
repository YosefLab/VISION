/*
   Initializes a zoomable scatter plot in the element "node"
   node = , for example, "#chart_div"
*/
function ColorScatter(node)
{
    this.node = $(node).get(0)
    this.width = $(node).width()
    this.height = $(node).height()
}

ColorScatter.prototype.clearData = function() {
}

/*
 * points: array of [x, y, color, name]
 * tree_points: array of [x, y]
 * tree_adj: array of [point1, points2] indices into tree_points
 *
 */
//ColorScatter.prototype.setData = function(points, isFactor, full_color_range, selected_cells, diverging_colormap)
ColorScatter.prototype.setData = function(object)
{
    var self = this
    var points = object['points']
    var isFactor = object['isFactor'] === undefined ? false : object['isFactor']
    var full_color_range = object['full_color_range'] === undefined ? false : object['full_color_range']
    var selected_cells = object['selected_cells'] === undefined ? [] : object['selected_cells']
    var diverging_colormap = object['diverging_colormap'] === undefined ? true : object['diverging_colormap']

    var circle_radius = this.pointsToRadius(points.length)

    var data = []; // Holds plotly traces

    var showlegend = false
    if(isFactor) {
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

        showlegend = true

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
            colorscale = [
                [0,   '#d8d8d8'],
                [0.5, '#395252'],
                [1,   '#000000'],
            ]
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
                showscale: true,
                colorbar: {
                    titleside: 'right', // default is 'top'
                    thickness: 15, // pixels
                    len: 0.5 //fraction of plot
                },
                colorscale: colorscale,
            },

        }

        data.push(trace1)

    }

    var shapes = []

    if ('tree_points' in object){
        _.forEach(object['tree_points'], tp => {
            var x = tp[0]
            var y = tp[1]
            var r = 5;
            shapes.push({
                type: 'circle',
                xref: 'x',
                yref: 'y',
                xsizemode: 'pixel',
                ysizemode: 'pixel',
                xanchor: x,
                yanchor: y,
                x0: -1*r,
                x1: r,
                y0: -1*r,
                y1: r,
                fillcolor: '#05ff65',
                opacity: 0.5,
                line: {
                    width: 1,
                },
            })
        })

        _.forEach(object['tree_adj'], tp => {
            var x0 = object['tree_points'][tp[0]][0]
            var y0 = object['tree_points'][tp[0]][1]
            var x1 = object['tree_points'][tp[1]][0]
            var y1 = object['tree_points'][tp[1]][1]
            shapes.push({
                type: 'line',
                xref: 'x',
                yref: 'y',
                x0: x0,
                x1: x1,
                y0: y0,
                y1: y1,
                line: {
                    color: "#05ff65",
                    width: 2,
                },
                opacity: 0.5,
            })
        })
    }


    var layout = {
        hovermode: 'closest',
        plot_bgcolor: '#eeeeee',
        showlegend: showlegend,
        legend: {
            xanchor: 'right',
            yanchor: 'right',
            x: 1,
            y: 1,
            bgcolor: 'rgba(255, 255, 255, .8)',
            bordercolor: 'rgba(0, 87, 82, .5)',
            borderwidth: 1,
        },
        margin: {
            l: 30,
            r: 90,
            t: 30,
            b: 30,
        },
        xaxis: {
            zeroline: false
        },
        yaxis: {
            zeroline: false
        },
        shapes: shapes,
    };

    var options = {
        'scrollZoom': true,
        'displaylogo': false,
        'displayModeBar': true,
        'modeBarButtonsToRemove': ['sendDataToCloud', 'hoverCompareCartesian', 'toggleSpikelines'],
    }

    Plotly.newPlot(this.node, data, layout, options)

    this.node.on('plotly_click', function(data){

        var point = data.points[0]
        var cellId = point.text

        var event = new CustomEvent('select-cells', {detail: [cellId]})
        window.dispatchEvent(event);
    })

    this.node.on('plotly_doubleclick', function(){
        // Unselect everything
        var selections = _.map(self.node.data, function(d){
            return _.range(d.x.length)
        });

        Plotly.restyle(self.node, {
            selectedpoints: selections,
        })

        var event = new CustomEvent('select-cells', {detail: []})
        window.dispatchEvent(event);
    })

    this.node.on('plotly_selected', function(eventData){
        var cellIds;
        if (eventData === undefined) {
            cellIds = [];
        } else {
            cellIds = _.map(eventData.points, d => {
                return d.data.text[d.pointNumber];
            });
        }
        var event = new CustomEvent('select-cells', {detail: cellIds})
        window.dispatchEvent(event);
    });

    this.node.on('plotly_legendclick', function(data){
        var selectedPoints = data.fullData[data.curveNumber].text

        // Update the selection in the plot
        var selections = _.map(data.fullData, function(d){
            if(d.index == data.curveNumber){
                return _.range(d.x.length)
            } else {
                return []
            }
        });

        Plotly.restyle(self.node, {
            selectedpoints: selections,
        })

        var event = new CustomEvent('select-cells', {detail: selectedPoints})
        window.dispatchEvent(event)
        return false;
    })

    this.node.on('plotly_legenddoubleclick', function(){
        return false;
    })
}

ColorScatter.prototype.pointsToRadius = function(n_points)
{
    // Pick a point size based on the number of scatter
    // plot points

    // Using a polynomial fit for this
    // Data taken with width: 1397, height: 790
    // Points      Low        High
    // 1052        8          12
    // 2884        7          9
    // 6857        5          7
    // 67171       2.5        3.5

    var scale_factor = Math.min(this.width, this.height)/790

    var x = 1/n_points
    var circle_radius =
        -1.3026e7*Math.pow(x, 2) +
        1.9784e4*Math.pow(x, 1) +
        2.9423*Math.pow(x, 0)

    circle_radius = Math.min(circle_radius, 12) // Looks silly if it's too big


    return circle_radius * scale_factor

}

/*
 * Automatically adjusts the scale if the points are too wide
 * or too narrow
 */
ColorScatter.prototype.autoZoom = function() {

};

ColorScatter.prototype.hover_cells = function(cell_ids) {

};

ColorScatter.prototype.resize = function()
{
    Plotly.Plots.resize(this.node);
};
