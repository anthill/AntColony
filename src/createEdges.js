'use strict'

var dt = require("delaunay-triangulate");

var points = require('./initializePoints.js').points;
var textMesh = require('./initializePoints.js').textMesh;
var textPointsId = require('./initializePoints.js').textPointsId;

var Edge = require('./edge.js');

// triangulate
var cells = dt(points.map(function(p){
    return [p.x, p.y]
}))

var edges = [];
var permutations = [[0,1], [0,2], [1,2]];

cells.forEach(function(cell){
   
    for (var i = 0; i < 3; ++i){  // for each point.id listed in current cell
        var pt = points[cell[i]];
        //console.log("-------- id:" + pt);

        for (var j = 1; j <= 2; ++j){ 

            var ptj = points[cell[( i + j ) % 3]]; // pick one of the other 2 points of the cell
            //console.log("--- other: " + ptj);
            var newEdge = undefined;

            // if pt already has nextEdges ...
            if (pt.nexts.length != 0) {
                //console.log("Trouvé 1");
                
                // ... get the points corresponding ...
                var tempPoints = pt.nexts.map(function(e){
                    return [e.pt1, e.pt2];
                })

                tempPoints = tempPoints.reduce(function(a, b){
                     return a.concat(b);
                });

                //console.log(tempPoints);

                // ... and check if ptj already is part of the existing nextEdges. If not, add the edge.
                if (tempPoints.indexOf(ptj) == -1){
                    //console.log("Trouvé 2");
                    //newEdge = createEdge(points[pt], points[ptj]);
                    newEdge = Edge.create(pt, ptj);
                    pt.nexts.push(newEdge);
                    edges.push(newEdge);
                    //console.log(edges.length);
                }
            }
            else {
                //console.log("Trouvé 3");
                newEdge = Edge.create(pt, ptj);
                pt.nexts.push(newEdge);
                edges.push(newEdge);
                //console.log(edges.length);
            }

            // add also the edge to the edge's other point's nextEdges
            if (newEdge != undefined){
                ptj.nexts.push(newEdge);
                //addToNextEdges(points[ptj], cell, temp, nextEdges); 
            }         
        }

        // add the textEdges to nextEdges map
        if (textMesh && textPointsId.indexOf(pt.id) != -1 && pt.id < (textPointsId.length - 1)) 
        {
            //console.log('verif');
            var textEdge = Edge.create(pt, points[pt.id + 1]);
            //var temp2 = [];
            edges.push(textEdge);
            pt.nexts.push(textEdge);
            //temp2.push(textEdge);
            //addToNextEdges(points[pt+1], undefined, temp2, nextEdges);
        }

        //console.log("longueur: " + points[pt].nexts.length);
        //addToNextEdges(points[pt], cell, nexts, nextEdges);
    }
    //console.log(edges.length);
})

module.exports = edges;