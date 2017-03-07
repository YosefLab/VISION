## Projection wrapper class


setMethod("initialize", signature(.Object = "ProjectionData"), 
          function(.Object, name, data) {
            
            .Object@name = name
            .Objection@data = data
            
            return(.Object)
          })