'use strict'

var dt = require("delaunay-triangulate");

var range = require('./utilities.js').range;

var points = require('./initializePoints.js').points;
var textMesh = require('./initializePoints.js').textMesh;
var citySet = require('./initializePoints.js').citySet;
var nbRandomPoints = require('./initializePoints.js').nbRandomPoints;
var forcedEdges = require('./initializePoints.js').forcedEdges;

var Edge = require('./edge.js');

// triangulate
var cells = dt(points.map(function(p){
    return [p.x, p.y]
}))

var edges = [];
var permutations = [[0,1], [0,2], [1,2]];

// force the edges of the text to be edges of the graph
if (textMesh) {
    range(0, points.length - nbRandomPoints).forEach(function(id){
        var directLink = forcedEdges[id];
        var textEdge = Edge.create(points[id], points[directLink]);
        edges.push(textEdge);
        points[id].nexts.push(textEdge);
    })
}


cells.forEach(function(cell){
   
    for (var i = 0; i < 3; ++i){  // for each point.id listed in current cell
        var pt = points[cell[i]];

        for (var j = 1; j < 3; ++j){ 

            var ptj = points[cell[( i + j ) % 3]]; // pick one of the other 2 points of the cell
            var newEdge = undefined;

            // if pt already has nextEdges ...
            if (pt.nexts.length != 0) {
                
                // ... get the points corresponding ...
                var tempPoints = pt.nexts.map(function(e){
                    return [e.pt1, e.pt2];
                }).reduce(function(a, b){
                     return a.concat(b);
                });

                // ... and check if ptj already is part of the existing nextEdges. If not, add the edge.
                if (tempPoints.indexOf(ptj) == -1){
                    newEdge = Edge.create(pt, ptj);
                    edges.push(newEdge);
                    pt.nexts.push(newEdge);
                }
            }
            else {
                newEdge = Edge.create(pt, ptj);
                edges.push(newEdge);
                pt.nexts.push(newEdge);
            }

            // add also the edge to the edge's other point's nextEdges
            if (newEdge != undefined){
                ptj.nexts.push(newEdge);
            }         
        }

        // add the textEdges to nextEdges map
        if (textMesh && (citySet.indexOf(pt) != -1)) {
            var textEdge = Edge.create(pt, points[pt.id + 1]);
            edges.push(textEdge);
            pt.nexts.push(textEdge);
        }

    }
})

module.exports = edges;