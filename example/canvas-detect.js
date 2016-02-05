'use strict';

/*
    Proudly copied from Modernizer
    https://github.com/Modernizr/Modernizr/blob/master/feature-detects/canvas.js
    MIT Licence
*/

module.exports = function(){
    var elem = document.createElement('canvas');
    return !!(elem.getContext && elem.getContext('2d'));
}