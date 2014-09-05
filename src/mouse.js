'use strict'

var mouse = {
    x: 0,
    y: 0,
    r: 0.03
};


window.addEventListener( 'mousemove', function(e){
    mouse.x = e.clientX / window.innerWidth;
    mouse.y = e.clientY / window.innerHeight;
    //var d = Math.sqrt(Math.pow(population[0].posX - mouse.x, 2) + Math.pow(population[0].posY - mouse.y, 2));
    
    // console.log('Souris: ' + mouse.x + '|' + mouse.y);
    // console.log('Ant: ' + population[0].posX + '|' + population[0].posY);
    // console.log('Distance: ' + d);
} );


module.exports = mouse;