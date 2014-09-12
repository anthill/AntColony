"use strict"

var points = require('./initializePoints.js').points;
var edges = require('./createEdges.js');
var population = require('./initializeAnts.js');
var shell = require('./rendering.js').shell;

var sqrt = Math.sqrt;
var pow = Math.pow;
var floor = Math.floor;
var random = Math.random;
var atan = Math.atan;

var norm = require('./utilities.js').norm;

// create points
var points = require('./initializePoints.js').points;

// create edges
var nbEdges = edges.length;

// initialize ants
var population = require('./initializeAnts.js');

// Rendering
var shell = require('./rendering.js').shell;
