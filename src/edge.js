'use strict';

//var Map = require("harmony-collections").Map;

var sqrt = Math.sqrt;
var pow = Math.pow;
var abs = Math.abs;
var atan = Math.atan;

//var norm = require('./utilities.js').norm;

var Vector = require('./vector.js');

// var nbEdges = 0;
// var nextEdges = new Map();


function Edge(ptA, ptB) {
    var distance = sqrt( pow(ptA.x - ptB.x, 2) + pow(ptA.y - ptB.y, 2) );
    //    direction = atan((ptB.y-ptA.y)/(ptB.x-ptA.x));

    // var direction;
    // if (ptA.x != ptB.x){
    //     direction = Math.atan((ptB.y-ptA.y)/(ptB.x-ptA.x));
    // } else {
    //     direction = Math.PI/2;
    // }

    // find line equation ax + by + c = 0
    var a = 1;
    var b = - (ptB.x - ptA.x) / (ptB.y - ptA.y);

    // orientate vector (a,b)
    if (b < 0){
        b = -b;
        a = -a;
    }

    // normalize vector (a,b)
    var n = new Vector(a, b);
    n.normalize();

    var c = - (a * ptA.x + b * ptA.y);

    // // calculate vector director
    var v = new Vector(ptB.x - ptA.x, ptB.y - ptA.y);
    
    v.normalize();


    this.id = undefined;
    this.pt1 = ptA;
    this.pt2 = ptB;
    this.direction = v;
    this.orthDirection = n; 
    this.distance = distance;
    this.pheromon = 1/distance;
    this.line = {
        a: a,
        b: b,
        c: c,
    }
}


// static methods
Edge.create = function(ptA, ptB) {
    var edge = new Edge(ptA, ptB);
    return edge;
}


// methods
Edge.prototype = {

    getOtherPoint: function(point) {
        if (point == this.pt1)
            return this.pt2;
        else if (point == this.pt2)
            return this.pt1;
        else
            console.log("Error");
    },

    calculateDistance: function(x, y) {
        var a = this.line.a,
            b = this.line.b,
            c = this.line.c;
        return abs(a * x + b * y + c) / Math.sqrt(Math.pow(a,2) + Math.pow(b,2));
    },

    /*calculateDistance: function(point) {
        var a = this.line.a,
            b = this.line.b,
            c = this.line.c;
        return abs(a * point.x + b * point.y + c) / Math.sqrt(Math.pow(a,2) + Math.pow(b,2));
    }*/
}
module.exports = Edge;