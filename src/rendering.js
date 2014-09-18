'use strict'

var random = Math.random;

module.exports = function(container){
    console.log('rendering', container);
    
    if(!container)
        throw new TypeError('Missing container');

    var points = require('./initializePoints.js').points;
    var citySet = require('./initializePoints.js').citySet;

    var edges = require('./createEdges.js');

    var population = require('./initializeAnts');
    var nbAnts = population.length;
        
    var canvas = document.createElement("canvas");
    canvas.width = window.innerWidth;
    canvas.height = window.innerHeight;
    container.appendChild(canvas);
    
    var context = canvas.getContext("2d");
    

    function tick() {
        var w = canvas.width;
        var h = canvas.height;
        var mouse = [lastMouseMoveEvent.clientX/w, lastMouseMoveEvent.clientY/h];
        context.setTransform(w, 0, 0, h, 0, 0);
        context.fillStyle = "#fff";
        context.fillRect(0,0,w,h);

        // edges
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

        // ants
        for(var i=0; i<population.length; ++i) {
            context.beginPath()
            var x = population[i].x + 0.005*random();
            var y = population[i].y + 0.005*random();
            // var x = population[i].x;
            // var y = population[i].y;
            context.fillStyle = "#1C1C1C"
            // context.arc(x, y, 0.002, 0, 2*Math.PI)
            context.fillRect(x, y, 0.0012, 0.0012);
            context.closePath();
            context.fill();
        }

        context.beginPath();
        context.fillStyle = "#AAAAAA";
        context.arc(mouse[0], mouse[1], 0.01, 0, 2*Math.PI);
        context.closePath();
        context.fill();
    };
    
    var lastMouseMoveEvent = {
        clientX: 0,
        clientY: 0
    };
    
    container.addEventListener('mousemove', function(e){
        lastMouseMoveEvent = e;
    });
    
    var paused = false;
    
    function togglePlayPause(){
        paused = !paused;
        if(!paused)
            animate();
    }
    
    container.addEventListener('click', togglePlayPause);
    
    function animate(){
        tick();
        
        if(!paused)
            requestAnimationFrame(animate);
    };
    animate();
    
    return {
        togglePlayPause: togglePlayPause,
        // should be a getter/setter, but IE8
        getAntCount: function(){
            throw 'TODO';
        },
        setAntCount: function(){
            throw 'TODO';
        }
    }
}