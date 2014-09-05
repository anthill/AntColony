'use strict'

// var nextEdges = require('./edge.js').nextEdges;
// var Edge = require('./edge.js').Edge;

function Point(x, y) {
    this.id = undefined;                
    this.x = x;
    this.y = y;
    this.nexts = [];
}

// Point.prototype.addToNextEdges = function(cell) {
    
//     if(nextEdges.has(this))
//     {
//         var existingNextEdges = nextEdges.get(this);
        
//         // for every already existing next edge
//         for (var i = 0; i < existingNextEdges.length; i++)
//         {
//             // check if current existing edge is part of current cell
//             if (cell != undefined && cell.indexOf(existingNextEdges[i].getOtherPoint(this)) == -1)
//                 // if not, i.e if current edge hasn't been already added, push current edge
//                 this.nexts.push(existingNextEdges[i]);
//         }   
//     }

//     // update nextEdges with new identified edges
//     nextEdges.set(this, this.nexts);

//     //console.log("NextEdges de: " + a.id);
//     //console.log(nextEdges.get(a));
// }

module.exports = Point;