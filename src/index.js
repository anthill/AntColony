"use strict"

var sqrt = Math.sqrt;
var pow = Math.pow;
var floor = Math.floor;
var random = Math.random;
var atan = Math.atan;

var norm = require('./utilities.js').norm;

// create points
var points = require('./initializePoints.js').points;

// create edges
var edges = require('./createEdges.js');
var nbEdges = edges.length;

// console.log(norm(edges[0].line.v));

// // edges verification
// for (var i = 0; i < points.length; i++) {
//     console.log('Longueur de ' + i + ': ' + points[i].nexts.length);
// }
// console.log('nb of edges: ' + nbEdges);

// initialize ants
var population = require('./initializeAnts.js');

// Rendering
var shell = require('./rendering.js').shell;