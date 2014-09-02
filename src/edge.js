'use strict';

var Map = require("harmony-collections").Map;

var sqrt = Math.sqrt;
var pow = Math.pow;
var atan = Math.atan;

// var nbEdges = 0;
// var nextEdges = new Map();


function Edge(ptA, ptB) {
    var distance = sqrt( pow(ptA.x - ptB.x, 2) + pow(ptA.y - ptB.y, 2) ),
        direction = atan((ptB.y-ptA.y)/(ptB.x-ptA.x));

    this.id = undefined;
    this.pt1 = ptA;
    this.pt2 = ptB;
    this.direction = direction; 
    this.distance = distance;
    this.pheromon = 1/distance;
}


// static methods
Edge.create = function(ptA, ptB) {
    var edge = new Edge(ptA, ptB);
    return edge;
}


// methods
Edge.prototype.getOtherPoint = function(point) {
    if (point == this.pt1)
        return this.pt2;
    else
        return this.pt1;
}

module.exports = Edge;