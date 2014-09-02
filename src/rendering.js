'use strict'

var shell = require("game-shell")();

var points = require('./initializePoints.js').points;
var citySet = require('./initializePoints.js').citySet;

var edges = require('./createEdges.js');

var population = require('./initializeAnts');
var nbAnts = population.length;

var canvas, context;

shell.on("init", function() {
    canvas = document.createElement("canvas");
    canvas.width = shell.width;
    canvas.height = shell.height;
    context = canvas.getContext("2d");
    shell.element.appendChild(canvas);
})

shell.on("resize", function(w, h) {
    canvas.width = w;
    canvas.height = h;
})

shell.on("render", function() {
    var w = canvas.width;
    var h = canvas.height;
    var mouse = [shell.mouseX/w, shell.mouseY/h];
    context.setTransform(w, 0, 0, h, 0, 0);
    context.fillStyle = "#fff";
    context.fillRect(0,0,w,h);

    // // edges
    // context.strokeStyle = "#000";
    // for(var i=0; i < edges.length; ++i) {
    //     context.lineWidth = 0.0001;
    //     var edge = edges[i];
    //     // if (edge.pheromon != 0){
    //     //     context.lineWidth = Math.min(0.00001 * edge.pheromon, 0.01);
    //     // } else {
    //     //     context.lineWidth = 0.00001;
    //     // }
    //     context.beginPath();
    //     context.moveTo(points[edge.pt1.id].x, points[edge.pt1.id].y);
    //     context.lineTo(points[edge.pt2.id].x, points[edge.pt2.id].y);
    //     context.stroke();
    // }

    // // vertices
    // for(var i=0; i<points.length; ++i) {
    //     context.beginPath()
    //     var point = points[i];
    //     if (citySet.has(point.id)) {
    //         context.fillStyle = "#0101DF";
    //         context.arc(point.x, point.y, 0.006, 0, 2*Math.PI);
    //     }
    //     else {
    //         context.fillStyle = "#000";
    //         context.arc(points[i].x, points[i].y, 0.003, 0, 2*Math.PI);
    //     }
    //     context.closePath();
    //     context.fill();
    // }

    // move ants
    for (i = 0; i < nbAnts; i++) {
        population[i].transit();
    }

    // pheromon evaporation
    for (i = 0; i < edges.length; i++) {
        if(edges[i].pheromon > 0){
            edges[i].pheromon -= 0.0001;
        }
    }


    for(var i=0; i<population.length; ++i) {
        context.beginPath()
        // var x = population[i].posX + 0.005*random();
        // var y = population[i].posY + 0.005*random();
        var x = population[i].posX;
        var y = population[i].posY;
        if (population[i].state === "pheromoning"){
            context.fillStyle = "#FF0000";
        }
        else {context.fillStyle = "#610B0B"}
        context.arc(x, y, 0.003, 0, 2*Math.PI)
        context.closePath()
        context.fill()
    }

    context.beginPath();
    context.fillStyle = "#AAAAAA"
    context.arc(mouse[0], mouse[1], 0.01, 0, 2*Math.PI);
    context.closePath();
    context.fill();
  
})

module.exports = {
	canvas: canvas,
	context: context,
	shell: shell
}