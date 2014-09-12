'use strict';

var random = Math.random;
var floor = Math.floor;

var sign = require('./utilities.js').sign;
var calculateDistance = require('./utilities.js').distance;
var norm = require('./utilities.js').norm;

var points = require('./initializePoints.js').points;
var citySet = require('./initializePoints.js').citySet;
var textPointsId = require('./initializePoints.js').textPointsId;
var possibleStartPointsId = require('./initializePoints.js').possibleStartPointsId;

var mouse = require('./mouse.js');


function Ant(point) {
    this.x = point.x;                
    this.y = point.y;
    this.velocity = 0.005;
    this.edge = undefined;
    //this.step = 0;
    this.state = "forage";
    this.edges = [];
    this.lastCity = undefined;
    this.origin = point;
    this.destination = undefined;
    this.orientation = undefined;
    this.prog = 0;
}
// forage: the ant wanders around without any pheromon deposition
// once it finds a city, it starts remembering the nodes it goes through
// when it finds another city, it computes the path length and adds pheromons one each edges
// proportionnaly to the shortestness of the path
// it resets the list of nodes and continues
// while foraging the ant choses the path with a pheromon preference


// static methods
Ant.generateRandStartPoint = function() {
    var randId = Math.floor(possibleStartPointsId.length * random());
    var randStartPoint = points[possibleStartPointsId[randId]];
    return randStartPoint;
}


// methods
Ant.prototype = {

    transit: function(){
        switch (this.state) {
        case "forage":
            var res = this.move();
            if (res.cityReached) {
                this.state = "pheromoning";
                this.lastCity = this.origin.id;
            };
            break;
        case "pheromoning":
            var res = this.move();
            if (res.edgeChanged) {
                this.edges.push(this.edge);
                // found a city
                if (res.cityReached && (this.origin.id != this.lastCity) ){
                    // compute the length of the path
                    var pathLength = this.edges.map(function(e){return e.distance}).reduce(function(a,b){return a + b});
                    var deltaPheromone = 1/pathLength;
                    this.edges.forEach(function(e){
                        var a = e.pt1, b = e.pt2, weight = 1;  
                        // increased dropped pheromons for textEdges
                        if (citySet.has(a.id) && citySet.has(b.id) && (Math.abs(a.id - b.id) == 1))
                        {
                            weight *= 10;
                        }
                        e.pheromon += (deltaPheromone * weight);
                    });

                    this.edges = [this.edge];
                    this.lastCity = this.origin.id;
                }
            }
          break;
        }

    },

    setDirection: function(){
        var possibleEdges = [];

        for (var i = 0; i < this.origin.nexts.length; i++)
        {
            possibleEdges[i] = this.origin.nexts[i];
        } 

        possibleEdges.splice(possibleEdges.indexOf(this.edge),1);

        // flip a coin and either take the smelliest path or a random one
        if (random() > 0.5){
            var smells = possibleEdges.map(function(e){return e.pheromon});
            var index = smells.indexOf(Math.max.apply(Math, smells));
            this.edge = possibleEdges[index];
        } 
        else
            this.edge = possibleEdges[floor(random()*possibleEdges.length)];

        // set the destination point, being edge.pt1 or edge.pt2
        this.destination = (this.origin == this.edge.pt1) ? this.edge.pt2 : this.edge.pt1;
        if (this.destination.x != this.origin.x)
            this.orientation = sign((this.destination.x-this.origin.x));
        else
            this.orientation = sign((this.destination.y-this.origin.y));
    },

    move: function(){
        var edgeChanged;
        var cityReached = false;

        // on edge
        if (this.prog < this.edge.distance){

            var delta = this.avoidObstacle();
            var distance = this.edge.calculateDistance(this.x, this.y);
            var deltaB = {x:0, y:0};

            // return to edge if too far
            if (distance > 0.003) {
                var line = this.edge.line;
                var sign = 1;

                if (this.y > -(line.a * this.x + line.c)/line.b)
                    sign *= -1;

                deltaB = {
                    x: sign * line.a,
                    y: sign * line.b
                }
            }

            this.x += (this.velocity*Math.cos(this.edge.direction)*this.orientation + delta.x*0.007 + deltaB.x*0.004);
            this.y += (this.velocity*Math.sin(this.edge.direction)*this.orientation + delta.y*0.007 + deltaB.y*0.004);

            this.prog = this.calculateProgression();
            // this.prog = calculateDistance(this, this.origin);
            console.log(this.prog/this.edge.distance);
            

            edgeChanged = false;

        // on vertex
        } else {
            this.step = 0;
            this.prog = 0;
            this.origin = this.destination;
            this.x = this.origin.x;
            this.y = this.origin.y;

            this.setDirection();

            cityReached = citySet.has(this.origin.id);
            edgeChanged = true;
        }
        return {cityReached: cityReached, edgeChanged: edgeChanged};
    },

    avoidObstacle: function(){
        // var distance = Math.sqrt(Math.pow(this.x - mouse.x, 2) + Math.pow(this.y - mouse.y, 2));
        var distance = calculateDistance(this, mouse);
    
        if (distance <= mouse.r)
        {
            return {
                x: (this.x - mouse.x)/distance,
                y: (this.y - mouse.y)/distance
                // x: (this.y - mouse.y)/distance,
                // y: - (this.x - mouse.x)/distance
            };

            this.direction
        }
        else
            return {x:0,y:0};
    },

    calculateProgression: function(){
        var v = {
            x: this.x - this.origin.x,
            y: this.y - this.origin.y
        }

        var theta = (v.x * this.edge.line.v.x + v.y * this.edge.line.v.y) / norm(v);
        var prog = norm(v) * Math.abs(theta)

        // returns length of projection on edge
        return prog;
    }

};

module.exports = Ant;
