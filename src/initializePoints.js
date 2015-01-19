'use strict';

var parse = require('parse-svg-path');

var range = require('./utilities.js').range;

var Point = require('./point.js');
var svgPath = require('./svgPath.js');

var random = Math.random;

var nbCity = 2;

var textMesh = true;

// Frame definition
var xInit = 0, yInit = 0;
var w = 1,
    h = 1;

var svgString = svgPath;

function svgToPoints(svgString) {
    var points = [];
    var edges = Object.create(null);

    var beginingPath;

    var X = 0;
    var Y = 0;
    var nbPoints = 0;
    var prevPoint;

    var commands = parse(svgString)
    for (var i=0; i<commands.length; i++){
        var command = commands[i];
        switch (command[0]) {
            case "m":
                X += command[1];
                Y += command[2];
                prevPoint = undefined;
                beginingPath = nbPoints;
                break;
            case "M":
                X = command[1];
                Y = command[2];
                prevPoint = undefined;
                beginingPath = nbPoints;
                break; 
            case "c":
                X = command[1];
                Y = command[2];
                points.push({id:nbPoints, x:X, y:Y});

                if (prevPoint != undefined) {
                    edges[prevPoint] = nbPoints;
                }
                prevPoint = nbPoints;
                nbPoints++;
                break; 
            case "l":
                X = command[1];
                Y = command[2];
                points.push({id:nbPoints, x:X, y:Y});

                if (prevPoint != undefined) {
                    edges[prevPoint] = nbPoints;
                }
                prevPoint = nbPoints;
                nbPoints++;
                break;
            case "z":
                edges[prevPoint] = nbPoints;
                beginingPath = undefined;
                prevPoint = undefined;
                break;    
        }
    }
    return {points : points, edges : edges};
}

// initialize points

module.exports = function(nbStartPoints, nbRandomPoints){
    var points = [];
    var forcedEdges;
    var citySet;

    if (textMesh){

        var myText = svgToPoints(svgString);
        points = myText.points;
        forcedEdges = myText.edges;
        citySet = range(0, points.length);

        var scaleX = 0.5;
        var scaleY = 0.5;
        var deltaX = 0.25;
        var deltaY = 0.2;

        // scale points to [0,1] + delta
        var maxX = Math.max.apply(Math, points.map(function(p){return p.x}));
        var minX = Math.min.apply(Math, points.map(function(p){return p.x}));
        var maxY = Math.max.apply(Math, points.map(function(p){return p.y}));
        var minY = Math.min.apply(Math, points.map(function(p){return p.y}));
        points = points.map(function(p){
            var x = scaleX * (p.x-minX)/(maxX-minX) + deltaX;
            var y = scaleY * (p.y-minY)/(maxY-minY) + deltaY;
            var newPoint = new Point(x, y);
            newPoint.id = p.id;

            return newPoint;
        });

        // only add random points
        var nbPoints = points.length;
        for(var i=0; i<nbRandomPoints; ++i) {

            var x = random();
            var y = random();

            var newPoint = new Point(x, y);
            newPoint.id = nbPoints;

            points.push(newPoint);

            nbPoints++;
        }

    } else {
        //add random points

        var nbPoints = 0;
        for(var i=0; i<nbRandomPoints; ++i) {

            var x = random();
            var y = random();

            var newPoint = new Point(x, y);
            newPoint.id = nbPoints;

            points.push(newPoint);
            
            nbPoints++;
        }

        citySet = range(0, nbCity);
        console.log(citySet);
    }


    // initialize start points
    var possibleStartPointsId = [];

    for (var i = 0; i < nbStartPoints; i++){
        possibleStartPointsId.push(Math.floor(nbPoints * random()));
    }
    

    return {
        textMesh: textMesh,
        points: points,
        citySet: citySet,
        possibleStartPointsId: possibleStartPointsId,
        nbRandomPoints: nbRandomPoints,
        forcedEdges: forcedEdges
    };
}
