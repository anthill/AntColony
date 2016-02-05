// 'use strict'; // This wouldn't allow js to be run in Safari :/

var isCanvasAvailable = require('./canvas-detect.js');
var terminal = require('./terminal.js');

var colonySection = document.querySelector('#colony');

var isDesktop = true;


// check if the device is mobile
var ua = window.navigator.userAgent;
if (ua.match(/Mobi/))
	isDesktop = false;

// don't run colony on mobiles
if (isDesktop){

	var antColony;

	if(isCanvasAvailable())
		antColony = terminal(colonySection);

	// // togglePlayPause
	// document.body.addEventListener('click', function(){
	// 	antColony.togglePlayPause();
	// });
}
