/*
   Initializes a zoomable scatter plot in the element "node"
   node = , for example, "#chart_div"
*/
function ColorScatter(node)
{
    this.node = $(node).get(0)
    this.firstPlot = true
    this.resize = _.debounce(this._resize, 300)
    this.title = ''
    this.initialZoom = null
    this.currentZoom = null
    this.n_points = null

    /*
     * This one is weird.  So when you double-click, plotly also
     * fires the click event, but fires it after the double-click on
     * some browsers (Chrome).  This causes the new point to be
     * selected where you double-clicked.  Apparently a known issue in Github.
     * To work around, we disable single-click events for a small
     * interval after a doubleclick using this variable
     */
    this.clickMask = false
}

/*
 * points: array of {'x', 'y', 'value', 'label', 'selected'}
 * tree_points: array of [x, y]
 * tree_adj: array of [point1, points2] indices into tree_points
 *
 */
ColorScatter.prototype.setData = function(object)
{
    var points = object['points']
    var isFactor = object['isFactor'] === undefined ? false : object['isFactor']
    var full_color_range = object['full_color_range'] === undefined ? false : object['full_color_range']
    var diverging_colormap = object['diverging_colormap'] === undefined ? true : object['diverging_colormap']
    var autoZoom = object['autozoom']

    this.title = object['title']
    this.n_points = points.length

    var x = _.map(points, p => p['x'])
    var xmin = _.min(x)
    var xmax = _.max(x)
    var xrange = xmax-xmin
    xmin = xmin - xrange * .05
    xmax = xmax + xrange * .05

    var y = _.map(points, p => p['y'])
    var ymin = _.min(y)
    var ymax = _.max(y)
    var yrange = ymax-ymin
    ymin = ymin - yrange * .05
    ymax = ymax + yrange * .05

    this.initialZoom = {
        'xmin': xmin,
        'xmax': xmax,
        'ymin': ymin,
        'ymax': ymax,
    }

    if(autoZoom || this.currentZoom === null) {

        this.currentZoom = {
            'xmin': xmin,
            'xmax': xmax,
            'ymin': ymin,
            'ymax': ymax,
        }
    }


    var circle_radius = this.pointsToRadius()

    var data = []; // Holds plotly traces

    var showlegend = false
    if(isFactor) {
        var c = _.map(points, p => p['value'])
        var unique = d3.set(c).values();

        // Need to check if any is selected.  If no, then all selected_points is null
        var anySelected = _(points).map('selected').reduce( (a, i) => a || i, false)

        _.forEach(unique, level => {
            var subset = _.filter(points, p => p['value'] == level)
            var x_sub = _.map(subset, p => p['x'])
            var y_sub = _.map(subset, p => p['y'])
            var id_sub = _.map(subset, p => p['label'])

            var selected_points
            if (anySelected) {
                selected_points = _(subset)
                    .map((p, i) => {return({selected: p.selected, index: i})})
                    .filter(p => p.selected)
                    .map(p => p.index)
                    .value()
            } else {
                selected_points = null
            }

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
                selectedpoints: selected_points,
                hoverinfo: 'text+name',
            }
            data.push(trace)
        });

        showlegend = true

    } else {

        var x = _.map(points, p => p['x'])
        var y = _.map(points, p => p['y'])
        var c = _.map(points, p => p['value'])
        var id = _.map(points, p => p['label'])

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

        var selected_points = _(points)
            .map((p, i) => {return({selected: p.selected, index: i})})
            .filter(p => p.selected)
            .map(p => p.index)
            .value()

        if(selected_points.length == 0)
            selected_points = null

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
            selectedpoints: selected_points,
            hoverinfo: 'text',
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

    var layout = this.getLayout()

    layout['shapes'] = shapes
    layout['showlegend'] = showlegend

    var options = this.getOptions()

    if (this.firstPlot || this.plotlyBug(data)) {
        Plotly.newPlot(this.node, data, layout, options)
        this.createListeners()
    } else {
        Plotly.react(this.node, data, layout, options)
    }

    this.firstPlot = false

}

ColorScatter.prototype.createListeners = function() {
    var self = this;

    this.node.on('plotly_click', function(data){

        if (self.clickMask) { return; }

        var point = data.points[0]
        var cellId = point.text

        var event = new CustomEvent('select-cells', {
            detail: { cells: [cellId] }
        })
        window.dispatchEvent(event);
    })

    this.node.on('plotly_doubleclick', function(){

        self.clickMask = true
        window.setTimeout(function() {self.clickMask = false; }, 500)

        var newLayout = {
            'xaxis.range[0]': self.initialZoom.xmin,
            'xaxis.range[1]': self.initialZoom.xmax,
            'yaxis.range[0]': self.initialZoom.ymin,
            'yaxis.range[1]': self.initialZoom.ymax,
        }

        self.relayout(newLayout)

        var event = new CustomEvent('select-cells', {
            detail: {cells: []}
        })
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
        var event = new CustomEvent('select-cells', {
            detail: {cells: cellIds}
        })
        window.dispatchEvent(event);
    });

    this.node.on('plotly_legendclick', function(data){
        var selectedPoints = data.fullData[data.curveNumber].text
        var name = data.fullData[data.curveNumber].name

        var event = new CustomEvent('select-cells', {
            detail: {cells: selectedPoints, name: name}
        })
        window.dispatchEvent(event)
        return false;
    })

    this.node.on('plotly_legenddoubleclick', function(){
        return false;
    })

    this.node.on('plotly_relayout', function(newLayout){

        if (_.isEmpty(newLayout)){ return; }

        var fireEvent = false
        var needsMarkerScale = false

        if (_.isEqual(newLayout, {autosize: true}))
            needsMarkerScale = true

        if ('xaxis.range[0]' in newLayout) {
            self.currentZoom.xmin = newLayout['xaxis.range[0]']
            fireEvent = true
            needsMarkerScale = true
        }

        if ('xaxis.range[1]' in newLayout) {
            self.currentZoom.xmax = newLayout['xaxis.range[1]']
            fireEvent = true
            needsMarkerScale = true
        }

        if ('yaxis.range[0]' in newLayout) {
            self.currentZoom.ymin = newLayout['yaxis.range[0]']
            fireEvent = true
            needsMarkerScale = true
        }

        if ('yaxis.range[1]' in newLayout) {
            self.currentZoom.ymax = newLayout['yaxis.range[1]']
            fireEvent = true
            needsMarkerScale = true
        }

        if (needsMarkerScale){
            var circle_radius = self.pointsToRadius()
            var dataUpdate = {
                'marker.size': circle_radius,
            }
            Plotly.restyle(self.node, dataUpdate)
        }

        if (fireEvent) {
            var event = new CustomEvent('scatter_relayout', {
                bubbles: true,
                detail: {newLayout: newLayout, origin: self},
            })
            self.node.dispatchEvent(event)
        }
    })
}

ColorScatter.prototype.getLayout = function() {

    var width = $(this.node).width()

    // Some options might depend on the plot size

    var titleOpt
    if(width < 800){
        titleOpt = {
            text: this.title,
            x: 0,
            xref: 'paper',
            xanchor: 'left',
            yanchor: 'top',
            font: {
                size: 14,
            }
        }

    } else {
        titleOpt = {
            text: this.title,
            x: 0,
            xref: 'paper',
            xanchor: 'left',
            yanchor: 'top',
            font: {
                size: 18,
            }
        }
    }

    var layout = {
        title: titleOpt,
        hovermode: 'closest',
        paper_bgcolor: 'rgba(255, 255, 255, 0)',
        plot_bgcolor: '#eeeeee',
        dragmode: 'pan',
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
            t: 50,
            b: 30,
        },
        xaxis: {
            zeroline: false,
            range: [this.currentZoom.xmin, this.currentZoom.xmax],
        },
        yaxis: {
            zeroline: false,
            range: [this.currentZoom.ymin, this.currentZoom.ymax],
        },
        modebar: {
            bgcolor: 'rgba(255, 255, 255, 0)',
        },
    }

    return layout
}

ColorScatter.prototype.getOptions = function()
{
    var width = $(this.node).width()
    var displayModeBar
    if(width < 800){
        displayModeBar = 'hover'
    } else {
        displayModeBar = true
    }

    var options = {
        'scrollZoom': true,
        'displaylogo': false,
        'displayModeBar': displayModeBar,
        'modeBarButtonsToRemove': ['sendDataToCloud', 'hoverCompareCartesian', 'toggleSpikelines'],
        'doubleClick': false,
    }

    return options

}

ColorScatter.prototype.updateSelection = function()
{
    var selected_cells = get_global_status('selected_cell')
    var selectedpoints;

    if(selected_cells.length <= 1){
        selectedpoints = [null]
    } else {
        var selected_cells_map = _.keyBy(selected_cells, x => x)
        selectedpoints = _(this.node.data)
            .map(trace => trace.text)
            .map(trace_ids => _(trace_ids)
                .map(x => x in selected_cells_map)
                .map((e, i) => {return({selected: e, index: i})}) // add indices
                .filter(x => x.selected)
                .map(x => x.index)
                .value()
            ).value()
    }

    // Unselect everything
    Plotly.restyle(this.node, {
        selectedpoints: selectedpoints,
    })
}

ColorScatter.prototype.pointsToRadius = function()
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

    var width = $(this.node).width()
    var height = $(this.node).height()
    var scale_factor = Math.min(width, height)/790

    var x = 1/this.n_points
    var circle_radius =
        -1.3026e7*Math.pow(x, 2) +
        1.9784e4*Math.pow(x, 1) +
        2.9423*Math.pow(x, 0)

    circle_radius = Math.min(circle_radius, 12) // Looks silly if it's too big

    circle_radius = circle_radius * scale_factor

    // Change markersize based on zoom factor
    var initialArea = (
        (this.initialZoom.xmax - this.initialZoom.xmin) *
        (this.initialZoom.ymax - this.initialZoom.ymin)
    )

    var currentArea = (
        (this.currentZoom.xmax - this.currentZoom.xmin) *
        (this.currentZoom.ymax - this.currentZoom.ymin)
    )

    circle_radius = circle_radius * Math.pow(initialArea/currentArea, 0.5)

    return circle_radius

}

/*
 * Automatically adjusts the scale if the points are too wide
 * or too narrow
 */
ColorScatter.prototype.autoZoom = function() {

};

ColorScatter.prototype.hover_cells = function(cell_ids) {

};

ColorScatter.prototype._resize = function()
{
    var layout = this.getLayout()
    Plotly.relayout(this.node, layout)

    Plotly.Plots.resize(this.node)
}

ColorScatter.prototype.relayout = function(newLayout)
{
    return Plotly.relayout(this.node, newLayout)
}

// Need to purge plots that aren't visible
// Or else there are too many WebGL contexts
ColorScatter.prototype.clear = function() {
    Plotly.purge(this.node)
    this.firstPlot = true
}

/*
 * Addresses Plotly Issue:
 *   - https://github.com/plotly/plotly.js/issues/3405
 */
ColorScatter.prototype.plotlyBug = function(newData) {
    var oldData = this.node.data

    var oldSizes = _.map(oldData, trace => trace.x.length)
    var newSizes = _.map(newData, trace => trace.x.length)

    var plotBug = false
    for(var i = 0; i < oldSizes.length; i++) {
        if ((oldSizes[i] > 100000) && (newSizes[i] <= 100000)) {
            plotBug = true
        }
    }

    return plotBug
}
