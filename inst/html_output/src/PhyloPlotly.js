// from https://github.com/jasondavies/newick.js/
function parseNewick(s) {
    var ancestors = [];
    var tree = {};
    var tokens = s.split(/\s*(;|\(|\)|,|:)\s*/);
    for (var i=0; i<tokens.length; i++) {
        var token = tokens[i];
        switch (token) {
            case '(': // new branchset
                var subtree = {};
                tree.branchset = [subtree];
                ancestors.push(tree);
                tree = subtree;
                break;
            case ',': // another branch
                var subtree = {};
                ancestors[ancestors.length-1].branchset.push(subtree);
                tree = subtree;
                break;
            case ')': // optional name next
                tree = ancestors.pop();
                break;
            case ':': // optional length next
                break;
            default:
                var x = tokens[i-1];
                if (x == ')' || x == '(' || x == ',') {
                    tree.name = token;
                } else if (x == ':') {
                    tree.length = parseFloat(token);
                }
        }
    }
    return tree;
}


/*
    Initialize a plotly phylo object with a newick sting and a div
    optionally, pass a dictionary of nodes -> values and a function to determine
    the colors of the the mapped values.
*/
function PhyloPlotly(newick, div, mapping, nodeColor) {
    this.newick = newick;
    this.div = div;
    this.tree = parseNewick(newick);
    this.mapping = mapping;
    this.nodeColor = nodeColor;
    this.collapsed = [];
    this.radial = true;
    this.depth = -1;
    
    // plotly values
    this.data = [];
    this.layout = {
        height:550,
        autoexpand: true,
        margin: {'l':10,'r':10,'t':10,'b':10},
        showlegend: false,
        yaxis: {'showgrid':false, 'zeroline':false, "showline":false, showticklabels:false},
        xaxis: {'showgrid':false, 'zeroline':false, "showline":false, showticklabels:false},
        hovermode:"closest"
    }
    
    this.yMin = 0;
    this.yMax = 1;
    this.d = 1
    
    this.collapseMethod = "mode";
    this.maxDepth = -1;
}


/*
  set the mapping
*/
PhyloPlotly.prototype.setMapping = function(mapping) {
  this.mapping = mapping;
}


/*
  set the collapseMethod for coloring
*/
PhyloPlotly.prototype.setCollapseMethod = function(collapseMethod) {
  this.collapseMethod = collapseMethod;
}

/*
  set the nodeColor function
*/
PhyloPlotly.prototype.setNodeColor = function(nodeColor) {
  this.nodeColor = nodeColor;
}


/*
  Init all the coordinates, plot for the first time
*/
PhyloPlotly.prototype.init = function() {
    this.updateCoords();
    this.generatePlotlyData();
    this.plotAndUpdateListeners();
}


/*
  Update the mode, linear or radial and redraw
*/
PhyloPlotly.prototype.updateMode = function(radial) {
    this.radial = radial;
    this.init();
}


/*
  Update depths and coordinates
*/
PhyloPlotly.prototype.updateCoords = function() {
    setUltrametricTreeDepths(this.tree, this.collapsed);
    this.setLinearCoords();
    this.setRadialCoords();
}


/*
  Populate the data array for plotly plotting
  item 0: horizontalLines, ie lines from root outwards
  item 1: vertical lines, ie lines perpendicular to root
  item 2: internal nodes and leaves
*/
PhyloPlotly.prototype.generatePlotlyData = function() {
    var xKey = "x";
    var yKey = "y";
    
    if (this.radial) {
        xKey = "u";
        yKey = "v";
    }
    
    // set the node coordinates
    var x = [];
    var y = [];
    var names = [];
    var colors = [];
    var nodeSizes = [];
    var nodeSymbols = [];
    
    var internalNodeSize = 5;
    var tipSize = 6.3;
    var collapsedNodeSize = 20;
    var self = this;
    
    function nodeCoords(tree) {
        names.push(tree["name"]);
        x.push(tree[xKey]);
        y.push(tree[yKey]);
        
        if (self.collapsed.includes(tree["name"])) {
            colors.push(self.nodeColor(self.collapseColor(tree)));
            nodeSizes.push(collapsedNodeSize);
            nodeSymbols.push("triangle-down");
            // black triangle for collapsed (color)
        } else if (!isTip(tree)) {
            colors.push("#a4a4a4");
            nodeSizes.push(internalNodeSize);
            nodeSymbols.push("circle");
            tree["branchset"].sort(function(a, b) {return a["depth"] - b["depth"]}).forEach(function(child) {
                nodeCoords(child);
            });
        } else {
            nodeSizes.push(tipSize);
            nodeSymbols.push("circle");
            colors.push(self.nodeColor(self.mapping[tree["name"]]));
        }
    }
    
    nodeCoords(this.tree);
    
    // set the lines
    
    var lines = this.linearLines();
    if (this.radial) {
        lines = this.radialLines();
    }
    
    var xh = lines[0];
    var yh = lines[1];
    var xv = lines[2];
    var yv = lines[3];

    var nodesTrace = {
        x: x,
        y: y,
        type: "scattergl",
        text: names,
        mode: "markers",
        marker: {
          color: colors,
          size: nodeSizes,
          symbol: nodeSymbols,
          opacity: 1,
          line: {width: 0}
        },
        hoverinfo: "text"
        
    }
      
    var horizontalLines = {
        x: xh,
        y: yh,
        mode: "lines",
        type:"scatter",
        line: {
            color: '#000000',
            width: 0.5,
            shape: "spline"
        },
        hoverinfo: 'skip'
    }
    
    var verticalLines = {
        x: xv,
        y: yv,
        mode: "lines",
        type:"scatter",
        line: {
            color: '#000000',
            width: 0.5,
            shape: "spline"
        },
        hoverinfo: 'skip'
    }
    this.data = [horizontalLines, verticalLines, nodesTrace]; 
}


/*
  Replot the tree and update the listeners.
*/
PhyloPlotly.prototype.plotAndUpdateListeners = function() {
    var self = this;
    $(self.div).off();
    Plotly.newPlot(self.div, self.data, self.layout);
    // add event listener for clicks
    
    // disable default right click on plotly div
    self.div.addEventListener('contextmenu', ev => {ev.preventDefault()});
    
    
    self.div.on('plotly_click', function(eventData) {
        var clickType = eventData.event.buttons;
        if (clickType === 2) {
          // right click
          // collapse or uncollapse all children
          var cell = eventData.points.map(x => x["text"])[0];
          var idx = self.collapsed.indexOf(cell);
          if (idx >= 0) {
            // remove cell (uncollapsing)
            self.collapsed.splice(idx, 1);
          } else {
            // add the cell (collapsing)
            self.collapsed.push(cell);
          }
          
          // update the coords then the data
          self.updateCoords();
          self.generatePlotlyData();
          // redraw the tree but with collapsed nodes as leaves
          Plotly.react(self.div, self.data, self.layout);
        } else {
          var expansion = selectionExpander(self.tree, eventData.points.map(x => x["text"]));
          var selectedNames = expansion.selectedNames;
          var cells = expansion.cells;
          
          var nodesToSelect = expansion.nodesToSelect;
          
          var selectedpointsIdx = [];
          var names = self.data[2].text;
          names.map(function (name, idx) {
              if (selectedNames.includes(name)) {
                  selectedpointsIdx.push(idx);
              }
          });
        
          self.data[2].selectedpoints =  selectedpointsIdx;
          
          Plotly.react(self.div, self.data, self.layout);
          
          set_global_status({"selected_cell":cells, "selection_type":"cells"});
        }
    });
    
    // add event listener for selection
    self.div.on('plotly_selected', function(eventData) {
        if (eventData === undefined) {
            set_global_status({"selected_cell":[], "selection_type":"none"});
            return;
        }
        // event data are points, internal and tips
        // for each tip use joining algo to find the internal nodes that lead to it
        // for each internal node find the tips and restyle so they're selected
        // restyle the lines so that they're highlight colored
        var expansion = selectionExpander(self.tree, eventData.points.map(x => x["text"]));
        var selectedNames = expansion.selectedNames;
        var cells = expansion.cells;
        
        var nodesToSelect = expansion.nodesToSelect;
        
        var selectedpointsIdx = [];
        var names = self.data[2].text;
        names.map(function (name, idx) {
            if (selectedNames.includes(name)) {
                selectedpointsIdx.push(idx);
            }
        });
      
        self.data[2].selectedpoints =  selectedpointsIdx;
        
        Plotly.react(self.div, self.data, self.layout);
        
        set_global_status({"selected_cell":cells, "selection_type":"cells"});
        
    })
}


// add the linear coords to the tree object
PhyloPlotly.prototype.setLinearCoords = function() {
    var maxDepth = this.tree["depth"]; // at root
    var leafIndex = 0;
    var y = [];
    function setLinearCoord(tree) {
        tree["x"] = maxDepth - tree["depth"]
        if (!isTip(tree)) {
            // internal node
            var childrenHeights = [];
            tree["branchset"].forEach(function(child) {
                setLinearCoord(child);
                childrenHeights.push(child["y"])
            });
            tree["y"] = (Math.max(...childrenHeights) + Math.min(...childrenHeights))/2;
        } else {
            // am a tip
            tree["y"] = leafIndex;
            leafIndex++;
            tree["x"] = tree["x"] + 1;
        }
        y.push(tree["y"]);
    }
    
    setLinearCoord(this.tree);
    this.yMin = Math.min(...y);
    this.yMax = Math.max(...y);
}


/* 
  add the radial coords to the tree object
*/
PhyloPlotly.prototype.setRadialCoords = function() {
    var self = this;
    function convertRadial(tree) {
        var radialCoords = linearCordToRadial([tree["x"], tree["y"]], self.yMax, self.yMin, self.d);
        var coords = polorCoordToCartestian(radialCoords);
        tree["r"] = radialCoords[0];
        tree["theta"] = radialCoords[1];
        tree["u"] = coords[0];
        tree["v"] = coords[1];
        if (self.collapsed.includes(tree["name"])) {
            // collapsed
        } else if (!isTip(tree)) {
            tree["branchset"].sort(function(a, b) {return a["depth"] - b["depth"]}).forEach(function(child) {
                convertRadial(child);
            });
        }
    }
    
    convertRadial(this.tree)
}


/*
  Generate the coordinates for the lines when in linear mode
*/
PhyloPlotly.prototype.linearLines = function() {
    var xh = [];
    var yh = [];
    var xv = [];
    var yv = [];
    var self = this;
    
    function lines(tree) {
        if (!isTip(tree)  && !self.collapsed.includes(tree["name"])) {
            var xStart = tree["x"];
            var y = tree["y"];
            
            var childrenXs = [];
            var childrenYs = [];
            tree["branchset"].sort(function(a, b) {
              return a["depth"] - b["depth"]
            }).forEach(function(child) {
                lines(child);
                childrenXs.push(child["x"]);
                childrenYs.push(child["y"]);
            });
            var xEnd = Math.min(...childrenXs);
            var hlineEnd = (xStart + xEnd)/2;
            xh = xh.concat([xStart, hlineEnd, null]) // first part of connection, x
            yh = yh.concat([y, y, null]) // first part of connection, y
            
            // draw the second parts horizontal x
            childrenXs.forEach(function(x) {
                xh = xh.concat([hlineEnd, x, null]);
            });
            
            // draw the second parts horizontal y
            childrenYs.forEach(function(y) {
                yh = yh.concat([y, y, null])
            })
            
            // draw the vertical part
            xv = xv.concat([hlineEnd, hlineEnd, null]);
            yv = yv.concat([Math.min(...childrenYs), Math.max(...childrenYs), null]);
        }
    }
    
    lines(self.tree);
    return([xh, yh, xv, yv]);
}


/*
  Generate the lines when in radial mode
*/
PhyloPlotly.prototype.radialLines = function() {
    var xh = [];
    var yh = [];
    var xv = [];
    var yv = [];
    var self = this;
    
    var rh = [];
    var thetah = [];
    var rv = []
    var thetav = [];
    function lines(tree) {
        if (!isTip(tree) && !self.collapsed.includes(tree["name"])) {
            var rStart = tree["r"]
            var theta = tree["theta"]
            
            var childrenRs = [];
            var childrenThetas = [];
            
            tree["branchset"].sort(function(a, b) {
              return a["theta"] - b["theta"]
            }).forEach(function(child) {
                lines(child);
                childrenRs.push(child["r"]);
                childrenThetas.push(child["theta"]);
            });
            var rEnd = Math.min(...childrenRs);
            var hlineEnd = (rStart + rEnd)/2;
            
            rh = rh.concat([rStart, hlineEnd, null]) // first part of connection, x
            thetah = thetah.concat([theta, theta, null]) // first part of connection, y
            
            
            // draw the second parts horizontal x
            childrenRs.forEach(function(r) {
                rh = rh.concat([hlineEnd, r, null]);
            });
            
            // draw the second parts horizontal y
            childrenThetas.forEach(function(theta) {
                thetah = thetah.concat([theta, theta, null])
            })
            
            // draw the vertical part as a spline
            var maxTheta = Math.max(...childrenThetas);
            var minTheta = Math.min(...childrenThetas);
            var thetaRange = _.range(minTheta, maxTheta, (maxTheta - minTheta)/25)
            childrenThetas.concat([theta]).concat(thetaRange).sort().forEach(function(thetaC) {
                rv.push(hlineEnd);
                thetav.push(thetaC);
            })
            
            rv.push(null);
            thetav.push(null);
            
        }
    }
    
    lines(self.tree);    
    
    var coords = _.unzip(_.zip(rh, thetah).map(polorCoordToCartestian));
    xh = coords[0];
    yh = coords[1];
    coords = _.unzip(_.zip(rv, thetav).map(polorCoordToCartestian));
    xv = coords[0];
    yv = coords[1];
    return([xh, yh, xv, yv]);
}


/*
  Update the selected CELLS from an arguments
*/
PhyloPlotly.prototype.updateSelection = function(cells) {
    var expansion = selectionExpander(this.tree, cells);
    var selectedNames = expansion.selectedNames;
    
    var names = this.data[2].text;
    var selectedpointsIdx = [];
    names.map(function (name, idx) {
        if (selectedNames.includes(name)) {
            selectedpointsIdx.push(idx);
        }
    });
    if (selectedpointsIdx.length > 0) {
        this.data[2].selectedpoints =  selectedpointsIdx;
    } else {
        this.data[2].selectedpoints =  null;
    }
    
    
    Plotly.react(this.div, this.data, this.layout);
}

/*
  Expand all Nodes
*/
PhyloPlotly.prototype.expandAll = function() {
  this.collapseDepth = 0;
  // reset the depths
  this.collapsed = [];
  setUltrametricTreeDepths(this.tree, this.collapsed);
  this.updateCoords();
  this.generatePlotlyData();
  Plotly.react(this.div, this.data, this.layout);
}


/*
  Collapse tree to certain depth and redraw
*/
PhyloPlotly.prototype.collapseToDepth = function(depth) {
  this.collapseDepth = depth;
  // reset the depths
  this.collapsed = [];
  setUltrametricTreeDepths(this.tree, this.collapsed);
  var max_depth = this.tree["depth"];
  this.maxDepth = max_depth - 1;
  depth = max_depth - depth;
  
  var self = this;
  
  // collapse to a specific detph
  function collapseAtDepth(tree) {
    if (tree["depth"] <= depth && !isTip(tree)) {
      self.collapsed.push(tree["name"]); 
    } else {
      if (!isTip(tree)) {
        tree["branchset"].forEach(function (x) {
          collapseAtDepth(x);
        });
      }
    }
  }
  
  collapseAtDepth(this.tree);
  this.updateCoords();
  self.generatePlotlyData();
  Plotly.react(this.div, this.data, this.layout);
}


/*
  function to color a collapsed node
*/
PhyloPlotly.prototype.collapseColor = function(tree) {
  var childValues = [];
  var self = this;
  function collectValues(tree) {
    if (!isTip(tree)) {
        tree["branchset"].forEach(function (x) {
          collectValues(x);
        });
      } else {
        childValues.push(self.mapping[tree["name"]])
      }
  }
  
  
  collectValues(tree);
  switch(this.collapseMethod) {
    case "arimean":
      return (ariMean(childValues));
    case "geomean":
      return (geoMean(childValues));
    case "median":
      return (median(childValues));
    case "mode":
      return (mode(childValues));
    default:
      console.log("defaulted to mode")
      return (mode(childValues));
  }
}


/*
  convert linear tree coords to radial r and theta
*/
function linearCordToRadial(coord, yMax, yMin, d) {
    var x = coord[0];
    var y = coord[1];
    
    if (x === null || y === null) {
        return [null, null]
    }
    
    var r = x;
    var theta = 2*Math.PI*(y-yMin+d)/(yMax-yMin+d);
    
    return [r, theta]
}


/*
  convert r and theta back to x and y
*/
function polorCoordToCartestian(polarCoord) {
    var r = polarCoord[0];
    var theta = polarCoord[1];
    
    if (r == null || theta == null) {
        return [null, null]
    }
    return [Math.cos(theta)*r, Math.sin(theta)*r]
}


/* 
  is this tree a tip?
*/
function isTip(tree) {
    return !(_.keys(tree).includes("branchset"));
}

/*
  set the modified ultrametric depths on a tree object
*/
function setUltrametricTreeDepths(tree, collapsed) {
    if (!isTip(tree) && !collapsed.includes(tree["name"])) {
        var my_depth = 0;
        tree["branchset"].forEach(function(child) {
            setUltrametricTreeDepths(child, collapsed);
            my_depth = Math.max(my_depth, child["depth"] + 1);
        });
        tree["depth"] = my_depth;
    } else {
        tree["depth"] = 1;
    }
}


/*
  Expand a selection of CELLS in a TREE to any internal nodes that are 
  fully satisfied by the leaves selected.
*/
function selectionExpander(tree, cells) {
    var selectedNames = cells;
    var cells = [];
    
    var nodesToSelect = [];
    function treeSelect(tree, parentSelected) {
        var updateSelected = parentSelected && !selectedNames.includes(tree["name"]);
        tree["selected"] = parentSelected || selectedNames.includes(tree["name"]);
        if (!isTip(tree)) {
            var childrenSelects = [];
            tree["branchset"].forEach(function(child) {
                treeSelect(child, tree["selected"]);
                childrenSelects.push(child["selected"]);
            });
            
            updateSelected = updateSelected || (!tree["selected"] && childrenSelects.every(x => x));
            tree["selected"] = tree["selected"] || childrenSelects.every(x => x);
        } else {
            if (tree["selected"]) {
                cells.push(tree["name"]);
            }
        }
        
        if (updateSelected) {
            selectedNames.push(tree["name"]);
            nodesToSelect.push(tree["name"]);
        }
    }
    
    treeSelect(tree, false);
    return({"cells":cells, "selectedNames":selectedNames, "nodesToSelect":nodesToSelect})
}


/* COLLAPSE COLORING METHODS
  -arithmetic mean
  -geometric mean
  -median
  -mode
*/

/*
  arithmetic mean of array
*/
function ariMean(arr) {
  return(arr.reduce(function(a, b) {return a+b}) / arr.length);
}

/*
  geometric mean of array
*/
function geoMean(arr) {
  return(Math.exp(arr.map(function(a) {return Math.log(a)}).reduce(function(a, b) {return a + b}) / arr.length));
}


/*
  median of an array
*/
function median(arr) {
  arr.sort();
  if (arr.length % 2 !== 0) {
    return(arr[arr.length / 2 - 0.5]);
  } else {
    return((arr[arr.length / 2] + arr[(arr.length / 2) - 1]) / 2);
  }
}


/*
  mode of an array
*/
function mode(arr) {
  var freq = {};
  arr.forEach(function(x) {
    if (x in freq) {
      freq[x] ++;
    } else {
      freq[x] = 1;
    }
  })
  return (Object.keys(freq).reduce((a, b) => freq[a] > freq[b] ? a : b));
}




/*
function plotDendro(newick, div, radial, mapping, nodeColor) {
    var tree = parseNewick(newick);
    
    function generatePlotlyData(collapsed) {
        // treat collapsed nodes as leaves (special leaves though)
        setUltrametricTreeDepths(tree, collapsed);
        generateLinearCoords(tree);
        
        var x = [];
        var y = [];
        var names = [];
        var colors = [];
        var nodeSizes = [];
        
        var internalNodeSize = 5;
        var tipSize = 6.3;
        
        function nodeCoords(tree) {
            names.push(tree["name"]);
            x.push(tree["x"]);
            y.push(tree["y"]);
            
            if (collapsed.includes(tree["name"])) {
              colors.push("#0000ff");
              nodeSizes.push(internalNodeSize);
              // blue triangle for collapsed (color)
            } else if (!isTip(tree)) {
                colors.push("#a4a4a4");
                nodeSizes.push(internalNodeSize);
                tree["branchset"].sort(function(a, b) {return a["depth"] - b["depth"]}).forEach(function(child) {
                    nodeCoords(child);
                });
            } else {
                nodeSizes.push(tipSize);
                colors.push(nodeColor(mapping[tree["name"]]));
            }
        }
        
        nodeCoords(tree);
        
        
        var yMin = Math.min(...y);
        var yMax = Math.max(...y);
        var d = 1;
        
        
        
        if (radial) {
            x = [];
            y = [];
            colors = [];
            nodeSizes = [];
            function convertRadial(tree) {
                var radialCoords = treeCordToRadial([tree["x"], tree["y"]]);
                var coords = polorCoordToCartestian(radialCoords);
                tree["r"] = radialCoords[0];
                tree["theta"] = radialCoords[1];
                tree["x"] = coords[0];
                tree["y"] = coords[1];
                x.push(coords[0]);
                y.push(coords[1]);
                if (collapsed.includes(tree["name"])) {
                    // collapsed
                    colors.push("#0000ff");
                    nodeSizes.push(internalNodeSize);
                } else if (!isTip(tree)) {
                    nodeSizes.push(internalNodeSize);
                    colors.push("#a4a4a4")
                    tree["branchset"].sort(function(a, b) {return a["depth"] - b["depth"]}).forEach(function(child) {
                        convertRadial(child);
                    });
                } else {
                    nodeSizes.push(tipSize);
                    colors.push(nodeColor(mapping[tree["name"]]));
                }
            }
            
            convertRadial(tree)
        }
        
        
        var xh = [];
        var yh = [];
        var xv = []
        var yv = [];
        
        if (radial) {
            var rh = [];
            var thetah = [];
            var rv = []
            var thetav = [];
            function lines(tree) {
                if (!isTip(tree) && !collapsed.includes(tree["name"])) {
                    var rStart = tree["r"]
                    var theta = tree["theta"]
                    
                    var childrenRs = [];
                    var childrenThetas = [];
                    
                    tree["branchset"].sort(function(a, b) {return a["theta"] - b["theta"]}).forEach(function(child) {
                        lines(child);
                        childrenRs.push(child["r"]);
                        childrenThetas.push(child["theta"]);
                    });
                    var rEnd = Math.min(...childrenRs);
                    var hlineEnd = (rStart + rEnd)/2;
                    
                    rh = rh.concat([rStart, hlineEnd, null]) // first part of connection, x
                    thetah = thetah.concat([theta, theta, null]) // first part of connection, y
                    
                    
                    // draw the second parts horizontal x
                    childrenRs.forEach(function(r) {
                        rh = rh.concat([hlineEnd, r, null]);
                    });
                    
                    // draw the second parts horizontal y
                    childrenThetas.forEach(function(theta) {
                        thetah = thetah.concat([theta, theta, null])
                    })
                    
                    // draw the vertical part as a spline
                    var maxTheta = Math.max(...childrenThetas);
                    var minTheta = Math.min(...childrenThetas);
                    var thetaRange = _.range(minTheta, maxTheta, (maxTheta - minTheta)/25)
                    childrenThetas.concat([theta]).concat(thetaRange).sort().forEach(function(thetaC) {
                        rv.push(hlineEnd);
                        thetav.push(thetaC);
                    })
                    
                    rv.push(null);
                    thetav.push(null);
                    
                }
            }
            
            lines(tree)    
            
            var coords = _.unzip(_.zip(rh, thetah).map(polorCoordToCartestian))
            xh = coords[0];
            yh = coords[1];
            coords = _.unzip(_.zip(rv, thetav).map(polorCoordToCartestian))
            xv = coords[0];
            yv = coords[1];
            
        } else {
            function lines(tree) {
                if (!isTip(tree)  && !collapsed.includes(tree["name"])) {
                    var xStart = tree["x"]
                    var y = tree["y"]
                    
                    var childrenXs = [];
                    var childrenYs = []
                    tree["branchset"].sort(function(a, b) {return a["depth"] - b["depth"]}).forEach(function(child) {
                        lines(child);
                        childrenXs.push(child["x"]);
                        childrenYs.push(child["y"]);
                    });
                    var xEnd = Math.min(...childrenXs);
                    var hlineEnd = (xStart + xEnd)/2;
                    xh = xh.concat([xStart, hlineEnd, null]) // first part of connection, x
                    yh = yh.concat([y, y, null]) // first part of connection, y
                    
                    // draw the second parts horizontal x
                    childrenXs.forEach(function(x) {
                        xh = xh.concat([hlineEnd, x, null]);
                    });
                    
                    // draw the second parts horizontal y
                    childrenYs.forEach(function(y) {
                        yh = yh.concat([y, y, null])
                    })
                    
                    // draw the vertical part
                    xv = xv.concat([hlineEnd, hlineEnd, null]);
                    yv = yv.concat([Math.min(...childrenYs), Math.max(...childrenYs), null]);
                }
            }
            lines(tree)
        }
        
        var nodesTrace = {
            x: x,
            y: y,
            type: "scattergl",
            text: names,
            mode: "markers",
            marker: {
              color: colors,
              size: nodeSizes,
              opacity: 1
            },
            hoverinfo: "text"
            
        }
        
        var horizontalLines = {
            x: xh,
            y: yh,
            mode: "lines",
            type:"scatter",
            line: {
              color: '#000000',
              width: 0.5,
              shape: "spline"
            },
            hoverinfo: 'skip'
        }
        
        var verticalLines = {
            x: xv,
            y: yv,
            mode: "lines",
            type:"scatter",
            line: {
              color: '#000000',
              width: 0.5,
              shape: "spline"
            },
            hoverinfo: 'skip'
        }
        return([horizontalLines, verticalLines, nodesTrace]);
    }
    
    
    
    
    
    var data = generatePlotlyData([]);
    // clear any existing event listeners
    
    return({"tree": tree, "data": data, "layout": layout})
}
*/
