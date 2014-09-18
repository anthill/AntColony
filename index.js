'use strict';

var initRendering = require('./src/rendering.js');

module.exports = function(containerElement){
    initRendering(containerElement);
    var points = require('./src/initializePoints.js');
    var edges = require('./src/createEdges.js');
};