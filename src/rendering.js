'use strict'

var antFunction = require('./ant.js');
var antsGroupFactory = require('./antsGroup');

var random = Math.random;

var RANDOMMVT = 0.003;
var ANTSIZE = 0.002;

module.exports = function(container, pointsMap, options){

	if(!container)
		throw new TypeError('Missing container');

	// Ants variables
	var edges = pointsMap.edges;
	var objPopulationMax = options.nbAnts;
	var objPopulationAdjusted = objPopulationMax;
	var pointsInfos = pointsMap.pointsInfos;
	var population = [];
	var nbAntsPerStep = 100;
	
	var Ant = antFunction(container, pointsInfos, options);
	var antsGroup = antsGroupFactory(Ant);

	// Animation variables
	var animID;
	var deltaTime;
	var FPSCount;
	var lastUpdate = performance.now();
	var FPSMonitor = document.querySelector('#FPS');
	var dTMonitor = document.querySelector('#dT');
	var warningMonitor = document.querySelector('#warning');
	var refreshTime = 0;
	var maxDeltaTime = 40;
	var FPSOverLimitCount = 0;
	var FPSUnderLimitCount = 0;

	// warning message disappears after 4 s
	warningMonitor.addEventListener("transitionend", function(){
		window.setTimeout(function(){
			warningMonitor.className = "invisible";
		}, 4000);
	})

	// Canvas
	var canvasList = document.getElementsByTagName("canvas");
	
	if (canvasList.length === 0){
		var canvas = document.createElement("canvas");
		var rect = container.getBoundingClientRect();
		canvas.width = rect.width;
		canvas.height = rect.height;
		canvas.style.backgroundColor = "rgb(255, 255, 255)"; 
		container.appendChild(canvas);
	}
	else{
		var canvas = canvasList[0];
		console.log('CANVAS');
	}
	var context = canvas.getContext("2d");
	context.clearRect ( 0 , 0 , canvas.width, canvas.height );
	

	function checkAntNumber(antNumber){
		if (antNumber < objPopulationAdjusted - 50){
			FPSMonitor.style.color = "green";
			population = antsGroup.create(population);
		}	
		else if (antNumber > objPopulationAdjusted){
			population = antsGroup.remove(population, antNumber - objPopulationAdjusted);
			FPSMonitor.style.color = "red";
			warningMonitor.className = "visible";
		}
		else{
			FPSMonitor.style.color = "white";
			// warningMonitor.className = "invisible";
		}
	}

	function displayFPS(dT){
		FPSCount = (1000/dT).toFixed(2);
		var t = dT.toFixed(2);
		FPSMonitor.textContent = 'FPS : ' + FPSCount;  
		dTMonitor.textContent = 'nbAnts : ' + population.length;
		// dTMonitor.innerText = 'dT : ' + t + 'ms';
	}

	function tick() {
		var now = performance.now();
		deltaTime = now - lastUpdate;
		lastUpdate = now;
		refreshTime += deltaTime/1000; // in seconds

		// console.log('nbAnts', population.length);

		checkAntNumber(population.length);

		// display FPS info every 0.3 s
		if (refreshTime > 0.3){
			displayFPS(deltaTime);
			refreshTime = 0; 
		}

		// remove ants when frame rate is too low
		if (FPSOverLimitCount === 10) {
			objPopulationAdjusted = objPopulationAdjusted * maxDeltaTime / deltaTime;
			FPSOverLimitCount = 0;
		}

		while (FPSUnderLimitCount > 50 && objPopulationAdjusted < objPopulationMax) {
			objPopulationAdjusted += 10;
		}

		// check duration of over/under framerate limit periods
		if (deltaTime > maxDeltaTime){
			FPSOverLimitCount++;
			FPSUnderLimitCount = 0;
		}
		else {
			FPSOverLimitCount = 0;
			FPSUnderLimitCount++;
		}

		// draw in canvas
		var w = canvas.width;
		var h = canvas.height;
		var mouse = [lastMouseMoveEvent.clientX/w, lastMouseMoveEvent.clientY/h];
		context.setTransform(w, 0, 0, h, 0, 0);
		context.fillStyle = "rgba(255, 255, 255, 0.6)";
		context.fillRect(0,0,w,h);

		// // edges
		// context.strokeStyle = "#000";
		// for(var i=0; i < edges.length; ++i) {
		//     context.lineWidth = 0.0001;
		//     var edge = edges[i];
		//     if (edge.pheromon != 0){
		//         context.lineWidth = Math.min(0.00001 * edge.pheromon, 0.01);
		//     } else {
		//         context.lineWidth = 0.00001;
		//     }
		//     context.beginPath();
		//     context.moveTo(pointsInfos.points[edge.pt1.id].x, pointsInfos.points[edge.pt1.id].y);
		//     context.lineTo(pointsInfos.points[edge.pt2.id].x, pointsInfos.points[edge.pt2.id].y);
		//     context.stroke();
		// }

		// // vertices
		// for(var i=0; i<pointsInfos.points.length; ++i) {
		//     context.beginPath()
		//     var point = pointsInfos.points[i];
		//     if (pointsInfos.citySet.indexOf(point.id) != -1){
		//         context.fillStyle = "#0101DF";
		//         context.arc(point.x, point.y, 0.006, 0, 2*Math.PI);
		//     }
		//     else {
		//         context.fillStyle = "#000";
		//         context.arc(pointsInfos.points[i].x, pointsInfos.points[i].y, 0.003, 0, 2*Math.PI);
		//     }
		//     context.closePath();
		//     context.fill();
		// }

		// move ants
		population.forEach(function(ant){
			ant.transit();
		});

		// pheromon evaporation
		edges.forEach(function(edge){
			if(edge.pheromon > 0){
				edge.pheromon -= 0.0001;
			}
		});

		// ants
		population.forEach(function(ant){
			context.beginPath()
			var x = ant.x + RANDOMMVT*random();
			var y = ant.y + RANDOMMVT*random();

			context.fillStyle = "black"
			context.fillRect(x, y, ANTSIZE, ANTSIZE);
			context.closePath();
			context.fill();
		})
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

	function reset(){
		population = [];
		edges = [];
		pointsInfos = [];

		cancelAnimationFrame(animID);
	}
	
	// container.addEventListener('click', togglePlayPause);

	function animate(){
		tick();
		
		if(!paused)
			animID = requestAnimationFrame(animate);
	}
	animate();

	function modifyAnts(opts){
		objPopulationMax = opts.nbAnts;
		objPopulationAdjusted = objPopulationMax;

		population.forEach(function(ant){
			ant.velocity = opts.velocity;
			ant.intelligence = opts.intelligence;
			ant.repSize = opts.repSize;
			ant.repSpeed = opts.repSpeed;
		});
	}
	
	return {
		togglePlayPause: togglePlayPause,
		reset: reset,
		// should be a getter/setter, but IE8
		getAntCount: function(){
			return population.length;
		},
		modifyAnts: modifyAnts
	}
}
