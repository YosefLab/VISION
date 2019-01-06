/*
A library for general purpose utility functions
*/


/*
 * Creates a detached copy of an SVG image
 * Copies all applied css styles directly
 * Used for exporting an SVG
 */
function svgCopy(svg)
{
    var svg_copy = svg.cloneNode(true);

    //Need to stick it in the DOM so basic styles are computed
    var div = document.createElement("div");
    div.style.display = "none";
    document.body.appendChild(div);
    div.appendChild(svg_copy);
    
    //Apply the differences here
    apply_styles(svg, svg_copy);
    return svg_copy;
}


/*
 * Recursively copy computed styles from one node to another
 * Repeats for all children
 * Assumes nodes between original and unstyled correspond directly
 */
function apply_styles(original_node, unstyled_node)
{
    var cssStyles = getComputedStyle(original_node);
    var unstyledStyles = getComputedStyle(unstyled_node);
    var computedStyleStr = unstyled_node.getAttribute("style");
    if(computedStyleStr === null){ computedStyleStr = "";}
    var i;
    for(i = 0; i < cssStyles.length;i++)
    {
        var key = cssStyles[i];
        var value = cssStyles.getPropertyValue(key);

        //For some reason, if height or width is set explicitly as an attribute, the computed style is "auto"
        //If we then add width:auto to the elements style string, we override the attribute
        //This prevents that from happening
        if(key === "height" && original_node.hasAttribute(key)){continue;}
        if(key === "width" && original_node.hasAttribute(key)){continue;}

        if (value !== unstyledStyles.getPropertyValue(key))
        {
            computedStyleStr += key + ":" + value + ";";
        }
    }

    unstyled_node.setAttribute("style", computedStyleStr);

    var original_children = original_node.children;
    var unstyled_children = unstyled_node.children;

    for(i=0; i < original_children.length; i++)
    {
        apply_styles(original_children[i], unstyled_children[i]);
    }
}


// Argsort an array in javascript
// Return the indices that would sort the array
Array.prototype.argSort = function()
{
    var out = new Array(this.length);
    for(var i = 0; i < out.length; i++) out[i] = i;
    var that = this;
    out.sort(function(a,b){return that[a] - that[b];});
    return out;
};

Array.prototype.argMin = function()
{
    var min = Infinity
    var min_i = -1
    for(var i = 0; i < this.length; i++){
        if(this[i] < min){
            min = this[i]
            min_i = i
        }
    }
    return min_i;
};

Array.prototype.argMax = function()
{
    var max = -Infinity
    var max_i = -1
    for(var i = 0; i < this.length; i++){
        if(this[i] > max){
            max = this[i]
            max_i = i
        }
    }
    return max_i;
};

var detect_browser_scrollbar_width = (function()
{
    var _width = -1; // cache the value so its only calculated once

    function inner()
    {
        if(_width > -1){
            return _width;
        }

        var n = $('<div/>')
            .css( {
                top: 0,
                left: $(window).scrollLeft()*-1,
                height: 1,
                width: 1,
                position: 'fixed',
                overflow: 'hidden'
            } )
            .append(
                $('<div/>')
                    .css( {
                        top: 1,
                        left: 1,
                        width: 150,
                        overflow: 'scroll',
                        position: 'absolute'
                    } )
                    .append(
                        $('<div/>')
                            .css( {
                                width: '100%',
                                height: 20
                            } )
                    )
            )
            .appendTo( 'body' );

        var child = n.children();

        var width = child[0].offsetWidth - child[0].clientWidth;

        n.remove()
        _width = width;

        return width
    }

    return inner;
}());


/*
 * Creates a stateful loading function
 *
 * Usage:
 *
 * var loadingFun = createLoadingFunction(my_div);
 *
 * loadingFun(true) // start loading after delay
 * loadingFun(false) // stop loading
 *
 */
function createLoadingFunction(node){
    var timer = -1;
    var loadingDiv;

    var loadFun = function(loadState){

        if(loadState === true){
            if(timer === -1){ // don't overwrite an existing timer and lose it
                timer = setTimeout(() => {
                    $(node).addClass('loading');
                    loadingDiv = document.createElement("div");
                    $(loadingDiv).height($(node).height());
                    $(loadingDiv).width($(node).width());
                    $(loadingDiv).offset($(node).offset());
                    $(loadingDiv).addClass('loadingSpinner');
                    var img = $('<img />', {
                        src: 'css/loading.svg',
                        alt: 'loading-spinner'
                    });
                    img.appendTo(loadingDiv);
                    $(node).parent().append(loadingDiv);
                }, 1000);
            }
        } else {
            if(timer !== -1){
                clearTimeout(timer);
                timer = -1;
            }
            $(node).removeClass('loading');
            $(loadingDiv).remove();
        }
    }

    return loadFun;
}
