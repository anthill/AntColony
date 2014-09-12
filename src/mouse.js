'use strict'

var Point = require('./point.js');
var edges = require('./createEdges.js');

var mouse = {
    x: 0,
    y: 0,
    r: 0.04
};


window.addEventListener( 'mousemove', function(e){
    mouse.x = e.clientX / window.innerWidth;
    mouse.y = e.clientY / window.innerHeight;
    //var d = Math.sqrt(Math.pow(population[0].posX - mouse.x, 2) + Math.pow(population[0].posY - mouse.y, 2));
    
    // console.log('Souris: ' + mouse.x + '|' + mouse.y);
    // console.log('Ant: ' + population[0].posX + '|' + population[0].posY);
    // console.log('Distance: ' + d);
});

window.addEventListener( 'click', function(e){
    var point = new Point(e.clientX / window.innerWidth, e.clientY / window.innerHeight);
    console.log(edges[0].calculateDistance(point));
});

module.exports = mouse;