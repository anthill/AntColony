(function e(t,n,r){function s(o,u){if(!n[o]){if(!t[o]){var a=typeof require=="function"&&require;if(!u&&a)return a(o,!0);if(i)return i(o,!0);var f=new Error("Cannot find module '"+o+"'");throw f.code="MODULE_NOT_FOUND",f}var l=n[o]={exports:{}};t[o][0].call(l.exports,function(e){var n=t[o][1][e];return s(n?n:e)},l,l.exports,e,t,n,r)}return n[o].exports}var i=typeof require=="function"&&require;for(var o=0;o<r.length;o++)s(r[o]);return s})({"/Users/Romain/Documents/Programmation/ants/AntColony/index.js":[function(require,module,exports){
'use strict';

var initRendering = require('./src/rendering.js');
var initializePoints = require('./src/initializePoints.js');
var createEdges = require('./src/createEdges.js');
// var initAnts = require('./src/initializeAnts');

module.exports = function init(containerElement, options){

	var render, pointsInfos, edges, population, pointsMap;


	function _init(containerElement, options){
		pointsInfos = initializePoints(options.nbStart, options.nbRand);
		edges = createEdges(pointsInfos.points);
		// population = options.nbAnts;
		// population = initAnts(containerElement, pointsInfos, options);
		pointsMap = {
			pointsInfos: pointsInfos,
			edges: edges
			// population: population
		};
		render = initRendering(containerElement, pointsMap, options);
	}

	_init(containerElement, options);

	return {
		togglePlayPause: function(){ render.togglePlayPause() },
		changeOptions: function(opts){
			render.modifyAnts(opts);
		},
		reset: function(opts){
			render.reset();

				// reset elements
			_init(containerElement, opts);
		}
	};
};
},{"./src/createEdges.js":"/Users/Romain/Documents/Programmation/ants/AntColony/src/createEdges.js","./src/initializePoints.js":"/Users/Romain/Documents/Programmation/ants/AntColony/src/initializePoints.js","./src/rendering.js":"/Users/Romain/Documents/Programmation/ants/AntColony/src/rendering.js"}],"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/ich.js":[function(require,module,exports){
"use strict"

//High level idea:
// 1. Use Clarkson's incremental construction to find convex hull
// 2. Point location in triangulation by jump and walk

module.exports = incrementalConvexHull

var orient = require("robust-orientation")
var compareCell = require("simplicial-complex").compareCells

function compareInt(a, b) {
  return a - b
}

function Simplex(vertices, adjacent, boundary) {
  this.vertices = vertices
  this.adjacent = adjacent
  this.boundary = boundary
  this.lastVisited = -1
}

Simplex.prototype.flip = function() {
  var t = this.vertices[0]
  this.vertices[0] = this.vertices[1]
  this.vertices[1] = t
  var u = this.adjacent[0]
  this.adjacent[0] = this.adjacent[1]
  this.adjacent[1] = u
}

function GlueFacet(vertices, cell, index) {
  this.vertices = vertices
  this.cell = cell
  this.index = index
}

function compareGlue(a, b) {
  return compareCell(a.vertices, b.vertices)
}

function bakeOrient(d) {
  var code = ["function orient(){var tuple=this.tuple;return test("]
  for(var i=0; i<=d; ++i) {
    if(i > 0) {
      code.push(",")
    }
    code.push("tuple[", i, "]")
  }
  code.push(")}return orient")
  var proc = new Function("test", code.join(""))
  var test = orient[d+1]
  if(!test) {
    test = orient
  }
  return proc(test)
}

var BAKED = []

function Triangulation(dimension, vertices, simplices) {
  this.dimension = dimension
  this.vertices = vertices
  this.simplices = simplices
  this.interior = simplices.filter(function(c) {
    return !c.boundary
  })

  this.tuple = new Array(dimension+1)
  for(var i=0; i<=dimension; ++i) {
    this.tuple[i] = this.vertices[i]
  }

  var o = BAKED[dimension]
  if(!o) {
    o = BAKED[dimension] = bakeOrient(dimension)
  }
  this.orient = o
}

var proto = Triangulation.prototype

//Degenerate situation where we are on boundary, but coplanar to face
proto.handleBoundaryDegeneracy = function(cell, point) {
  var d = this.dimension
  var n = this.vertices.length - 1
  var tuple = this.tuple
  var verts = this.vertices

  //Dumb solution: Just do dfs from boundary cell until we find any peak, or terminate
  var toVisit = [ cell ]
  cell.lastVisited = -n
  while(toVisit.length > 0) {
    cell = toVisit.pop()
    var cellVerts = cell.vertices
    var cellAdj = cell.adjacent
    for(var i=0; i<=d; ++i) {
      var neighbor = cellAdj[i]
      if(!neighbor.boundary || neighbor.lastVisited <= -n) {
        continue
      }
      var nv = neighbor.vertices
      for(var j=0; j<=d; ++j) {
        var vv = nv[j]
        if(vv < 0) {
          tuple[j] = point
        } else {
          tuple[j] = verts[vv]
        }
      }
      var o = this.orient()
      if(o > 0) {
        return neighbor
      }
      neighbor.lastVisited = -n
      if(o === 0) {
        toVisit.push(neighbor)
      }
    }
  }
  return null
}

proto.walk = function(point, random) {
  //Alias local properties
  var n = this.vertices.length - 1
  var d = this.dimension
  var verts = this.vertices
  var tuple = this.tuple

  //Compute initial jump cell
  var initIndex = random ? (this.interior.length * Math.random())|0 : (this.interior.length-1)
  var cell = this.interior[ initIndex ]

  //Start walking
outerLoop:
  while(!cell.boundary) {
    var cellVerts = cell.vertices
    var cellAdj = cell.adjacent

    for(var i=0; i<=d; ++i) {
      tuple[i] = verts[cellVerts[i]]
    }
    cell.lastVisited = n

    //Find farthest adjacent cell
    for(var i=0; i<=d; ++i) {
      var neighbor = cellAdj[i]
      if(neighbor.lastVisited >= n) {
        continue
      }
      var prev = tuple[i]
      tuple[i] = point
      var o = this.orient()
      tuple[i] = prev
      if(o < 0) {
        cell = neighbor
        continue outerLoop
      } else {
        if(!neighbor.boundary) {
          neighbor.lastVisited = n
        } else {
          neighbor.lastVisited = -n
        }
      }
    }
    return
  }

  return cell
}

proto.addPeaks = function(point, cell) {
  var n = this.vertices.length - 1
  var d = this.dimension
  var verts = this.vertices
  var tuple = this.tuple
  var interior = this.interior
  var simplices = this.simplices

  //Walking finished at boundary, time to add peaks
  var tovisit = [ cell ]

  //Stretch initial boundary cell into a peak
  cell.lastVisited = n
  cell.vertices[cell.vertices.indexOf(-1)] = n
  cell.boundary = false
  interior.push(cell)

  //Record a list of all new boundaries created by added peaks so we can glue them together when we are all done
  var glueFacets = []

  //Do a traversal of the boundary walking outward from starting peak
  while(tovisit.length > 0) {
    //Pop off peak and walk over adjacent cells
    var cell = tovisit.pop()
    var cellVerts = cell.vertices
    var cellAdj = cell.adjacent
    var indexOfN = cellVerts.indexOf(n)
    if(indexOfN < 0) {
      continue
    }

    for(var i=0; i<=d; ++i) {
      if(i === indexOfN) {
        continue
      }

      //For each boundary neighbor of the cell
      var neighbor = cellAdj[i]
      if(!neighbor.boundary || neighbor.lastVisited >= n) {
        continue
      }

      var nv = neighbor.vertices

      //Test if neighbor is a peak
      if(neighbor.lastVisited !== -n) {      
        //Compute orientation of p relative to each boundary peak
        var indexOfNeg1 = 0
        for(var j=0; j<=d; ++j) {
          if(nv[j] < 0) {
            indexOfNeg1 = j
            tuple[j] = point
          } else {
            tuple[j] = verts[nv[j]]
          }
        }
        var o = this.orient()

        //Test if neighbor cell is also a peak
        if(o > 0) {
          nv[indexOfNeg1] = n
          neighbor.boundary = false
          interior.push(neighbor)
          tovisit.push(neighbor)
          neighbor.lastVisited = n
          continue
        } else {
          neighbor.lastVisited = -n
        }
      }

      var na = neighbor.adjacent

      //Otherwise, replace neighbor with new face
      var vverts = cellVerts.slice()
      var vadj = cellAdj.slice()
      var ncell = new Simplex(vverts, vadj, true)
      simplices.push(ncell)

      //Connect to neighbor
      var opposite = na.indexOf(cell)
      if(opposite < 0) {
        continue
      }
      na[opposite] = ncell
      vadj[indexOfN] = neighbor

      //Connect to cell
      vverts[i] = -1
      vadj[i] = cell
      cellAdj[i] = ncell

      //Flip facet
      ncell.flip()

      //Add to glue list
      for(var j=0; j<=d; ++j) {
        var uu = vverts[j]
        if(uu < 0 || uu === n) {
          continue
        }
        var nface = new Array(d-1)
        var nptr = 0
        for(var k=0; k<=d; ++k) {
          var vv = vverts[k]
          if(vv < 0 || k === j) {
            continue
          }
          nface[nptr++] = vv
        }
        glueFacets.push(new GlueFacet(nface, ncell, j))
      }
    }
  }

  //Glue boundary facets together
  glueFacets.sort(compareGlue)

  for(var i=0; i+1<glueFacets.length; i+=2) {
    var a = glueFacets[i]
    var b = glueFacets[i+1]
    var ai = a.index
    var bi = b.index
    if(ai < 0 || bi < 0) {
      continue
    }
    a.cell.adjacent[a.index] = b.cell
    b.cell.adjacent[b.index] = a.cell
  }
}

proto.insert = function(point, random) {
  //Add point
  var verts = this.vertices
  verts.push(point)

  var cell = this.walk(point, random)
  if(!cell) {
    return
  }

  //Alias local properties
  var d = this.dimension
  var tuple = this.tuple

  //Degenerate case: If point is coplanar to cell, then walk until we find a non-degenerate boundary
  for(var i=0; i<=d; ++i) {
    var vv = cell.vertices[i]
    if(vv < 0) {
      tuple[i] = point
    } else {
      tuple[i] = verts[vv]
    }
  }
  var o = this.orient(tuple)
  if(o < 0) {
    return
  } else if(o === 0) {
    cell = this.handleBoundaryDegeneracy(cell, point)
    if(!cell) {
      return
    }
  }

  //Add peaks
  this.addPeaks(point, cell)
}

//Extract all boundary cells
proto.boundary = function() {
  var d = this.dimension
  var boundary = []
  var cells = this.simplices
  var nc = cells.length
  for(var i=0; i<nc; ++i) {
    var c = cells[i]
    if(c.boundary) {
      var bcell = new Array(d)
      var cv = c.vertices
      var ptr = 0
      var parity = 0
      for(var j=0; j<=d; ++j) {
        if(cv[j] >= 0) {
          bcell[ptr++] = cv[j]
        } else {
          parity = j&1
        }
      }
      if(parity === (d&1)) {
        var t = bcell[0]
        bcell[0] = bcell[1]
        bcell[1] = t
      }
      boundary.push(bcell)
    }
  }
  return boundary
}

function incrementalConvexHull(points, randomSearch) {
  var n = points.length
  if(n === 0) {
    throw new Error("Must have at least d+1 points")
  }
  var d = points[0].length
  if(n <= d) {
    throw new Error("Must input at least d+1 points")
  }

  //FIXME: This could be degenerate, but need to select d+1 non-coplanar points to bootstrap process
  var initialSimplex = points.slice(0, d+1)

  //Make sure initial simplex is positively oriented
  var o = orient.apply(void 0, initialSimplex)
  if(o === 0) {
    throw new Error("Input not in general position")
  }
  var initialCoords = new Array(d+1)
  for(var i=0; i<=d; ++i) {
    initialCoords[i] = i
  }
  if(o < 0) {
    initialCoords[0] = 1
    initialCoords[1] = 0
  }

  //Create initial topological index, glue pointers together (kind of messy)
  var initialCell = new Simplex(initialCoords, new Array(d+1), false)
  var boundary = initialCell.adjacent
  var list = new Array(d+2)
  for(var i=0; i<=d; ++i) {
    var verts = initialCoords.slice()
    for(var j=0; j<=d; ++j) {
      if(j === i) {
        verts[j] = -1
      }
    }
    var t = verts[0]
    verts[0] = verts[1]
    verts[1] = t
    var cell = new Simplex(verts, new Array(d+1), true)
    boundary[i] = cell
    list[i] = cell
  }
  list[d+1] = initialCell
  for(var i=0; i<=d; ++i) {
    var verts = boundary[i].vertices
    var adj = boundary[i].adjacent
    for(var j=0; j<=d; ++j) {
      var v = verts[j]
      if(v < 0) {
        adj[j] = initialCell
        continue
      }
      for(var k=0; k<=d; ++k) {
        if(boundary[k].vertices.indexOf(v) < 0) {
          adj[j] = boundary[k]
        }
      }
    }
  }

  //Initialize triangles
  var triangles = new Triangulation(d, initialSimplex, list)

  //Insert remaining points
  var useRandom = !!randomSearch
  for(var i=d+1; i<n; ++i) {
    triangles.insert(points[i], useRandom)
  }
  
  //Extract boundary cells
  return triangles.boundary()
}
},{"robust-orientation":"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/orientation.js","simplicial-complex":"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/simplicial-complex/topology.js"}],"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-scale/node_modules/two-sum/two-sum.js":[function(require,module,exports){
"use strict"

module.exports = fastTwoSum

function fastTwoSum(a, b, result) {
	var x = a + b
	var bv = x - a
	var av = x - bv
	var br = b - bv
	var ar = a - av
	if(result) {
		result[0] = ar + br
		result[1] = x
		return result
	}
	return [ar+br, x]
}
},{}],"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-scale/robust-scale.js":[function(require,module,exports){
"use strict"

var twoProduct = require("two-product")
var twoSum = require("two-sum")

module.exports = scaleLinearExpansion

function scaleLinearExpansion(e, scale) {
  var n = e.length
  if(n === 1) {
    var ts = twoProduct(e[0], scale)
    if(ts[0]) {
      return ts
    }
    return [ ts[1] ]
  }
  var g = new Array(2 * n)
  var q = [0.1, 0.1]
  var t = [0.1, 0.1]
  var count = 0
  twoProduct(e[0], scale, q)
  if(q[0]) {
    g[count++] = q[0]
  }
  for(var i=1; i<n; ++i) {
    twoProduct(e[i], scale, t)
    var pq = q[1]
    twoSum(pq, t[0], q)
    if(q[0]) {
      g[count++] = q[0]
    }
    var a = t[1]
    var b = q[1]
    var x = a + b
    var bv = x - a
    var y = b - bv
    q[1] = x
    if(y) {
      g[count++] = y
    }
  }
  if(q[1]) {
    g[count++] = q[1]
  }
  if(count === 0) {
    g[count++] = 0.0
  }
  g.length = count
  return g
}
},{"two-product":"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/two-product/two-product.js","two-sum":"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-scale/node_modules/two-sum/two-sum.js"}],"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-subtract/robust-diff.js":[function(require,module,exports){
"use strict"

module.exports = robustSubtract

//Easy case: Add two scalars
function scalarScalar(a, b) {
  var x = a + b
  var bv = x - a
  var av = x - bv
  var br = b - bv
  var ar = a - av
  var y = ar + br
  if(y) {
    return [y, x]
  }
  return [x]
}

function robustSubtract(e, f) {
  var ne = e.length|0
  var nf = f.length|0
  if(ne === 1 && nf === 1) {
    return scalarScalar(e[0], -f[0])
  }
  var n = ne + nf
  var g = new Array(n)
  var count = 0
  var eptr = 0
  var fptr = 0
  var abs = Math.abs
  var ei = e[eptr]
  var ea = abs(ei)
  var fi = -f[fptr]
  var fa = abs(fi)
  var a, b
  if(ea < fa) {
    b = ei
    eptr += 1
    if(eptr < ne) {
      ei = e[eptr]
      ea = abs(ei)
    }
  } else {
    b = fi
    fptr += 1
    if(fptr < nf) {
      fi = -f[fptr]
      fa = abs(fi)
    }
  }
  if((eptr < ne && ea < fa) || (fptr >= nf)) {
    a = ei
    eptr += 1
    if(eptr < ne) {
      ei = e[eptr]
      ea = abs(ei)
    }
  } else {
    a = fi
    fptr += 1
    if(fptr < nf) {
      fi = -f[fptr]
      fa = abs(fi)
    }
  }
  var x = a + b
  var bv = x - a
  var y = b - bv
  var q0 = y
  var q1 = x
  var _x, _bv, _av, _br, _ar
  while(eptr < ne && fptr < nf) {
    if(ea < fa) {
      a = ei
      eptr += 1
      if(eptr < ne) {
        ei = e[eptr]
        ea = abs(ei)
      }
    } else {
      a = fi
      fptr += 1
      if(fptr < nf) {
        fi = -f[fptr]
        fa = abs(fi)
      }
    }
    b = q0
    x = a + b
    bv = x - a
    y = b - bv
    if(y) {
      g[count++] = y
    }
    _x = q1 + x
    _bv = _x - q1
    _av = _x - _bv
    _br = x - _bv
    _ar = q1 - _av
    q0 = _ar + _br
    q1 = _x
  }
  while(eptr < ne) {
    a = ei
    b = q0
    x = a + b
    bv = x - a
    y = b - bv
    if(y) {
      g[count++] = y
    }
    _x = q1 + x
    _bv = _x - q1
    _av = _x - _bv
    _br = x - _bv
    _ar = q1 - _av
    q0 = _ar + _br
    q1 = _x
    eptr += 1
    if(eptr < ne) {
      ei = e[eptr]
    }
  }
  while(fptr < nf) {
    a = fi
    b = q0
    x = a + b
    bv = x - a
    y = b - bv
    if(y) {
      g[count++] = y
    } 
    _x = q1 + x
    _bv = _x - q1
    _av = _x - _bv
    _br = x - _bv
    _ar = q1 - _av
    q0 = _ar + _br
    q1 = _x
    fptr += 1
    if(fptr < nf) {
      fi = -f[fptr]
    }
  }
  if(q0) {
    g[count++] = q0
  }
  if(q1) {
    g[count++] = q1
  }
  if(!count) {
    g[count++] = 0.0  
  }
  g.length = count
  return g
}
},{}],"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-sum/robust-sum.js":[function(require,module,exports){
"use strict"

module.exports = linearExpansionSum

//Easy case: Add two scalars
function scalarScalar(a, b) {
  var x = a + b
  var bv = x - a
  var av = x - bv
  var br = b - bv
  var ar = a - av
  var y = ar + br
  if(y) {
    return [y, x]
  }
  return [x]
}

function linearExpansionSum(e, f) {
  var ne = e.length|0
  var nf = f.length|0
  if(ne === 1 && nf === 1) {
    return scalarScalar(e[0], f[0])
  }
  var n = ne + nf
  var g = new Array(n)
  var count = 0
  var eptr = 0
  var fptr = 0
  var abs = Math.abs
  var ei = e[eptr]
  var ea = abs(ei)
  var fi = f[fptr]
  var fa = abs(fi)
  var a, b
  if(ea < fa) {
    b = ei
    eptr += 1
    if(eptr < ne) {
      ei = e[eptr]
      ea = abs(ei)
    }
  } else {
    b = fi
    fptr += 1
    if(fptr < nf) {
      fi = f[fptr]
      fa = abs(fi)
    }
  }
  if((eptr < ne && ea < fa) || (fptr >= nf)) {
    a = ei
    eptr += 1
    if(eptr < ne) {
      ei = e[eptr]
      ea = abs(ei)
    }
  } else {
    a = fi
    fptr += 1
    if(fptr < nf) {
      fi = f[fptr]
      fa = abs(fi)
    }
  }
  var x = a + b
  var bv = x - a
  var y = b - bv
  var q0 = y
  var q1 = x
  var _x, _bv, _av, _br, _ar
  while(eptr < ne && fptr < nf) {
    if(ea < fa) {
      a = ei
      eptr += 1
      if(eptr < ne) {
        ei = e[eptr]
        ea = abs(ei)
      }
    } else {
      a = fi
      fptr += 1
      if(fptr < nf) {
        fi = f[fptr]
        fa = abs(fi)
      }
    }
    b = q0
    x = a + b
    bv = x - a
    y = b - bv
    if(y) {
      g[count++] = y
    }
    _x = q1 + x
    _bv = _x - q1
    _av = _x - _bv
    _br = x - _bv
    _ar = q1 - _av
    q0 = _ar + _br
    q1 = _x
  }
  while(eptr < ne) {
    a = ei
    b = q0
    x = a + b
    bv = x - a
    y = b - bv
    if(y) {
      g[count++] = y
    }
    _x = q1 + x
    _bv = _x - q1
    _av = _x - _bv
    _br = x - _bv
    _ar = q1 - _av
    q0 = _ar + _br
    q1 = _x
    eptr += 1
    if(eptr < ne) {
      ei = e[eptr]
    }
  }
  while(fptr < nf) {
    a = fi
    b = q0
    x = a + b
    bv = x - a
    y = b - bv
    if(y) {
      g[count++] = y
    } 
    _x = q1 + x
    _bv = _x - q1
    _av = _x - _bv
    _br = x - _bv
    _ar = q1 - _av
    q0 = _ar + _br
    q1 = _x
    fptr += 1
    if(fptr < nf) {
      fi = f[fptr]
    }
  }
  if(q0) {
    g[count++] = q0
  }
  if(q1) {
    g[count++] = q1
  }
  if(!count) {
    g[count++] = 0.0  
  }
  g.length = count
  return g
}
},{}],"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/two-product/two-product.js":[function(require,module,exports){
"use strict"

module.exports = twoProduct

var SPLITTER = +(Math.pow(2, 27) + 1.0)

function twoProduct(a, b, result) {
  var x = a * b

  var c = SPLITTER * a
  var abig = c - a
  var ahi = c - abig
  var alo = a - ahi

  var d = SPLITTER * b
  var bbig = d - b
  var bhi = d - bbig
  var blo = b - bhi

  var err1 = x - (ahi * bhi)
  var err2 = err1 - (alo * bhi)
  var err3 = err2 - (ahi * blo)

  var y = alo * blo - err3

  if(result) {
    result[0] = y
    result[1] = x
    return result
  }

  return [ y, x ]
}
},{}],"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/orientation.js":[function(require,module,exports){
"use strict"

var twoProduct = require("two-product")
var robustSum = require("robust-sum")
var robustScale = require("robust-scale")
var robustSubtract = require("robust-subtract")

var NUM_EXPAND = 5

var EPSILON     = 1.1102230246251565e-16
var ERRBOUND3   = (3.0 + 16.0 * EPSILON) * EPSILON
var ERRBOUND4   = (7.0 + 56.0 * EPSILON) * EPSILON

function cofactor(m, c) {
  var result = new Array(m.length-1)
  for(var i=1; i<m.length; ++i) {
    var r = result[i-1] = new Array(m.length-1)
    for(var j=0,k=0; j<m.length; ++j) {
      if(j === c) {
        continue
      }
      r[k++] = m[i][j]
    }
  }
  return result
}

function matrix(n) {
  var result = new Array(n)
  for(var i=0; i<n; ++i) {
    result[i] = new Array(n)
    for(var j=0; j<n; ++j) {
      result[i][j] = ["m", j, "[", (n-i-1), "]"].join("")
    }
  }
  return result
}

function sign(n) {
  if(n & 1) {
    return "-"
  }
  return ""
}

function generateSum(expr) {
  if(expr.length === 1) {
    return expr[0]
  } else if(expr.length === 2) {
    return ["sum(", expr[0], ",", expr[1], ")"].join("")
  } else {
    var m = expr.length>>1
    return ["sum(", generateSum(expr.slice(0, m)), ",", generateSum(expr.slice(m)), ")"].join("")
  }
}

function determinant(m) {
  if(m.length === 2) {
    return [["sum(prod(", m[0][0], ",", m[1][1], "),prod(-", m[0][1], ",", m[1][0], "))"].join("")]
  } else {
    var expr = []
    for(var i=0; i<m.length; ++i) {
      expr.push(["scale(", generateSum(determinant(cofactor(m, i))), ",", sign(i), m[0][i], ")"].join(""))
    }
    return expr
  }
}

function orientation(n) {
  var pos = []
  var neg = []
  var m = matrix(n)
  var args = []
  for(var i=0; i<n; ++i) {
    if((i&1)===0) {
      pos.push.apply(pos, determinant(cofactor(m, i)))
    } else {
      neg.push.apply(neg, determinant(cofactor(m, i)))
    }
    args.push("m" + i)
  }
  var posExpr = generateSum(pos)
  var negExpr = generateSum(neg)
  var funcName = "orientation" + n + "Exact"
  var code = ["function ", funcName, "(", args.join(), "){var p=", posExpr, ",n=", negExpr, ",d=sub(p,n);\
return d[d.length-1];};return ", funcName].join("")
  var proc = new Function("sum", "prod", "scale", "sub", code)
  return proc(robustSum, twoProduct, robustScale, robustSubtract)
}

var orientation3Exact = orientation(3)
var orientation4Exact = orientation(4)

var CACHED = [
  function orientation0() { return 0 },
  function orientation1() { return 0 },
  function orientation2(a, b) { 
    return b[0] - a[0]
  },
  function orientation3(a, b, c) {
    var l = (a[1] - c[1]) * (b[0] - c[0])
    var r = (a[0] - c[0]) * (b[1] - c[1])
    var det = l - r
    var s
    if(l > 0) {
      if(r <= 0) {
        return det
      } else {
        s = l + r
      }
    } else if(l < 0) {
      if(r >= 0) {
        return det
      } else {
        s = -(l + r)
      }
    } else {
      return det
    }
    var tol = ERRBOUND3 * s
    if(det >= tol || det <= -tol) {
      return det
    }
    return orientation3Exact(a, b, c)
  },
  function orientation4(a,b,c,d) {
    var adx = a[0] - d[0]
    var bdx = b[0] - d[0]
    var cdx = c[0] - d[0]
    var ady = a[1] - d[1]
    var bdy = b[1] - d[1]
    var cdy = c[1] - d[1]
    var adz = a[2] - d[2]
    var bdz = b[2] - d[2]
    var cdz = c[2] - d[2]
    var bdxcdy = bdx * cdy
    var cdxbdy = cdx * bdy
    var cdxady = cdx * ady
    var adxcdy = adx * cdy
    var adxbdy = adx * bdy
    var bdxady = bdx * ady
    var det = adz * (bdxcdy - cdxbdy) 
            + bdz * (cdxady - adxcdy)
            + cdz * (adxbdy - bdxady)
    var permanent = (Math.abs(bdxcdy) + Math.abs(cdxbdy)) * Math.abs(adz)
                  + (Math.abs(cdxady) + Math.abs(adxcdy)) * Math.abs(bdz)
                  + (Math.abs(adxbdy) + Math.abs(bdxady)) * Math.abs(cdz)
    var tol = ERRBOUND4 * permanent
    if ((det > tol) || (-det > tol)) {
      return det
    }
    return orientation4Exact(a,b,c,d)
  }
]

function slowOrient(args) {
  var proc = CACHED[args.length]
  if(!proc) {
    proc = CACHED[args.length] = orientation(args.length)
  }
  return proc.apply(undefined, args)
}

function generateOrientationProc() {
  while(CACHED.length <= NUM_EXPAND) {
    CACHED.push(orientation(CACHED.length))
  }
  var args = []
  var procArgs = ["slow"]
  for(var i=0; i<=NUM_EXPAND; ++i) {
    args.push("a" + i)
    procArgs.push("o" + i)
  }
  var code = [
    "function getOrientation(", args.join(), "){switch(arguments.length){case 0:case 1:return 0;"
  ]
  for(var i=2; i<=NUM_EXPAND; ++i) {
    code.push("case ", i, ":return o", i, "(", args.slice(0, i).join(), ");")
  }
  code.push("}var s=new Array(arguments.length);for(var i=0;i<arguments.length;++i){s[i]=arguments[i]};return slow(s);}return getOrientation")
  procArgs.push(code.join(""))

  var proc = Function.apply(undefined, procArgs)
  module.exports = proc.apply(undefined, [slowOrient].concat(CACHED))
  for(var i=0; i<=NUM_EXPAND; ++i) {
    module.exports[i] = CACHED[i]
  }
}

generateOrientationProc()
},{"robust-scale":"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-scale/robust-scale.js","robust-subtract":"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-subtract/robust-diff.js","robust-sum":"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-sum/robust-sum.js","two-product":"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/two-product/two-product.js"}],"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/simplicial-complex/node_modules/bit-twiddle/twiddle.js":[function(require,module,exports){
/**
 * Bit twiddling hacks for JavaScript.
 *
 * Author: Mikola Lysenko
 *
 * Ported from Stanford bit twiddling hack library:
 *    http://graphics.stanford.edu/~seander/bithacks.html
 */

"use strict"; "use restrict";

//Number of bits in an integer
var INT_BITS = 32;

//Constants
exports.INT_BITS  = INT_BITS;
exports.INT_MAX   =  0x7fffffff;
exports.INT_MIN   = -1<<(INT_BITS-1);

//Returns -1, 0, +1 depending on sign of x
exports.sign = function(v) {
  return (v > 0) - (v < 0);
}

//Computes absolute value of integer
exports.abs = function(v) {
  var mask = v >> (INT_BITS-1);
  return (v ^ mask) - mask;
}

//Computes minimum of integers x and y
exports.min = function(x, y) {
  return y ^ ((x ^ y) & -(x < y));
}

//Computes maximum of integers x and y
exports.max = function(x, y) {
  return x ^ ((x ^ y) & -(x < y));
}

//Checks if a number is a power of two
exports.isPow2 = function(v) {
  return !(v & (v-1)) && (!!v);
}

//Computes log base 2 of v
exports.log2 = function(v) {
  var r, shift;
  r =     (v > 0xFFFF) << 4; v >>>= r;
  shift = (v > 0xFF  ) << 3; v >>>= shift; r |= shift;
  shift = (v > 0xF   ) << 2; v >>>= shift; r |= shift;
  shift = (v > 0x3   ) << 1; v >>>= shift; r |= shift;
  return r | (v >> 1);
}

//Computes log base 10 of v
exports.log10 = function(v) {
  return  (v >= 1000000000) ? 9 : (v >= 100000000) ? 8 : (v >= 10000000) ? 7 :
          (v >= 1000000) ? 6 : (v >= 100000) ? 5 : (v >= 10000) ? 4 :
          (v >= 1000) ? 3 : (v >= 100) ? 2 : (v >= 10) ? 1 : 0;
}

//Counts number of bits
exports.popCount = function(v) {
  v = v - ((v >>> 1) & 0x55555555);
  v = (v & 0x33333333) + ((v >>> 2) & 0x33333333);
  return ((v + (v >>> 4) & 0xF0F0F0F) * 0x1010101) >>> 24;
}

//Counts number of trailing zeros
function countTrailingZeros(v) {
  var c = 32;
  v &= -v;
  if (v) c--;
  if (v & 0x0000FFFF) c -= 16;
  if (v & 0x00FF00FF) c -= 8;
  if (v & 0x0F0F0F0F) c -= 4;
  if (v & 0x33333333) c -= 2;
  if (v & 0x55555555) c -= 1;
  return c;
}
exports.countTrailingZeros = countTrailingZeros;

//Rounds to next power of 2
exports.nextPow2 = function(v) {
  v += v === 0;
  --v;
  v |= v >>> 1;
  v |= v >>> 2;
  v |= v >>> 4;
  v |= v >>> 8;
  v |= v >>> 16;
  return v + 1;
}

//Rounds down to previous power of 2
exports.prevPow2 = function(v) {
  v |= v >>> 1;
  v |= v >>> 2;
  v |= v >>> 4;
  v |= v >>> 8;
  v |= v >>> 16;
  return v - (v>>>1);
}

//Computes parity of word
exports.parity = function(v) {
  v ^= v >>> 16;
  v ^= v >>> 8;
  v ^= v >>> 4;
  v &= 0xf;
  return (0x6996 >>> v) & 1;
}

var REVERSE_TABLE = new Array(256);

(function(tab) {
  for(var i=0; i<256; ++i) {
    var v = i, r = i, s = 7;
    for (v >>>= 1; v; v >>>= 1) {
      r <<= 1;
      r |= v & 1;
      --s;
    }
    tab[i] = (r << s) & 0xff;
  }
})(REVERSE_TABLE);

//Reverse bits in a 32 bit word
exports.reverse = function(v) {
  return  (REVERSE_TABLE[ v         & 0xff] << 24) |
          (REVERSE_TABLE[(v >>> 8)  & 0xff] << 16) |
          (REVERSE_TABLE[(v >>> 16) & 0xff] << 8)  |
           REVERSE_TABLE[(v >>> 24) & 0xff];
}

//Interleave bits of 2 coordinates with 16 bits.  Useful for fast quadtree codes
exports.interleave2 = function(x, y) {
  x &= 0xFFFF;
  x = (x | (x << 8)) & 0x00FF00FF;
  x = (x | (x << 4)) & 0x0F0F0F0F;
  x = (x | (x << 2)) & 0x33333333;
  x = (x | (x << 1)) & 0x55555555;

  y &= 0xFFFF;
  y = (y | (y << 8)) & 0x00FF00FF;
  y = (y | (y << 4)) & 0x0F0F0F0F;
  y = (y | (y << 2)) & 0x33333333;
  y = (y | (y << 1)) & 0x55555555;

  return x | (y << 1);
}

//Extracts the nth interleaved component
exports.deinterleave2 = function(v, n) {
  v = (v >>> n) & 0x55555555;
  v = (v | (v >>> 1))  & 0x33333333;
  v = (v | (v >>> 2))  & 0x0F0F0F0F;
  v = (v | (v >>> 4))  & 0x00FF00FF;
  v = (v | (v >>> 16)) & 0x000FFFF;
  return (v << 16) >> 16;
}


//Interleave bits of 3 coordinates, each with 10 bits.  Useful for fast octree codes
exports.interleave3 = function(x, y, z) {
  x &= 0x3FF;
  x  = (x | (x<<16)) & 4278190335;
  x  = (x | (x<<8))  & 251719695;
  x  = (x | (x<<4))  & 3272356035;
  x  = (x | (x<<2))  & 1227133513;

  y &= 0x3FF;
  y  = (y | (y<<16)) & 4278190335;
  y  = (y | (y<<8))  & 251719695;
  y  = (y | (y<<4))  & 3272356035;
  y  = (y | (y<<2))  & 1227133513;
  x |= (y << 1);
  
  z &= 0x3FF;
  z  = (z | (z<<16)) & 4278190335;
  z  = (z | (z<<8))  & 251719695;
  z  = (z | (z<<4))  & 3272356035;
  z  = (z | (z<<2))  & 1227133513;
  
  return x | (z << 2);
}

//Extracts nth interleaved component of a 3-tuple
exports.deinterleave3 = function(v, n) {
  v = (v >>> n)       & 1227133513;
  v = (v | (v>>>2))   & 3272356035;
  v = (v | (v>>>4))   & 251719695;
  v = (v | (v>>>8))   & 4278190335;
  v = (v | (v>>>16))  & 0x3FF;
  return (v<<22)>>22;
}

//Computes next combination in colexicographic order (this is mistakenly called nextPermutation on the bit twiddling hacks page)
exports.nextCombination = function(v) {
  var t = v | (v - 1);
  return (t + 1) | (((~t & -~t) - 1) >>> (countTrailingZeros(v) + 1));
}


},{}],"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/simplicial-complex/node_modules/union-find/index.js":[function(require,module,exports){
"use strict"; "use restrict";

module.exports = UnionFind;

function UnionFind(count) {
  this.roots = new Array(count);
  this.ranks = new Array(count);
  
  for(var i=0; i<count; ++i) {
    this.roots[i] = i;
    this.ranks[i] = 0;
  }
}

var proto = UnionFind.prototype

Object.defineProperty(proto, "length", {
  "get": function() {
    return this.roots.length
  }
})

proto.makeSet = function() {
  var n = this.roots.length;
  this.roots.push(n);
  this.ranks.push(0);
  return n;
}

proto.find = function(x) {
  var roots = this.roots;
  while(roots[x] !== x) {
    var y = roots[x];
    roots[x] = roots[y];
    x = y;
  }
  return x;
}

proto.link = function(x, y) {
  var xr = this.find(x)
    , yr = this.find(y);
  if(xr === yr) {
    return;
  }
  var ranks = this.ranks
    , roots = this.roots
    , xd    = ranks[xr]
    , yd    = ranks[yr];
  if(xd < yd) {
    roots[xr] = yr;
  } else if(yd < xd) {
    roots[yr] = xr;
  } else {
    roots[yr] = xr;
    ++ranks[xr];
  }
}
},{}],"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/simplicial-complex/topology.js":[function(require,module,exports){
"use strict"; "use restrict";

var bits      = require("bit-twiddle")
  , UnionFind = require("union-find")

//Returns the dimension of a cell complex
function dimension(cells) {
  var d = 0
    , max = Math.max
  for(var i=0, il=cells.length; i<il; ++i) {
    d = max(d, cells[i].length)
  }
  return d-1
}
exports.dimension = dimension

//Counts the number of vertices in faces
function countVertices(cells) {
  var vc = -1
    , max = Math.max
  for(var i=0, il=cells.length; i<il; ++i) {
    var c = cells[i]
    for(var j=0, jl=c.length; j<jl; ++j) {
      vc = max(vc, c[j])
    }
  }
  return vc+1
}
exports.countVertices = countVertices

//Returns a deep copy of cells
function cloneCells(cells) {
  var ncells = new Array(cells.length)
  for(var i=0, il=cells.length; i<il; ++i) {
    ncells[i] = cells[i].slice(0)
  }
  return ncells
}
exports.cloneCells = cloneCells

//Ranks a pair of cells up to permutation
function compareCells(a, b) {
  var n = a.length
    , t = a.length - b.length
    , min = Math.min
  if(t) {
    return t
  }
  switch(n) {
    case 0:
      return 0;
    case 1:
      return a[0] - b[0];
    case 2:
      var d = a[0]+a[1]-b[0]-b[1]
      if(d) {
        return d
      }
      return min(a[0],a[1]) - min(b[0],b[1])
    case 3:
      var l1 = a[0]+a[1]
        , m1 = b[0]+b[1]
      d = l1+a[2] - (m1+b[2])
      if(d) {
        return d
      }
      var l0 = min(a[0], a[1])
        , m0 = min(b[0], b[1])
        , d  = min(l0, a[2]) - min(m0, b[2])
      if(d) {
        return d
      }
      return min(l0+a[2], l1) - min(m0+b[2], m1)
    
    //TODO: Maybe optimize n=4 as well?
    
    default:
      var as = a.slice(0)
      as.sort()
      var bs = b.slice(0)
      bs.sort()
      for(var i=0; i<n; ++i) {
        t = as[i] - bs[i]
        if(t) {
          return t
        }
      }
      return 0
  }
}
exports.compareCells = compareCells

function compareZipped(a, b) {
  return compareCells(a[0], b[0])
}

//Puts a cell complex into normal order for the purposes of findCell queries
function normalize(cells, attr) {
  if(attr) {
    var len = cells.length
    var zipped = new Array(len)
    for(var i=0; i<len; ++i) {
      zipped[i] = [cells[i], attr[i]]
    }
    zipped.sort(compareZipped)
    for(var i=0; i<len; ++i) {
      cells[i] = zipped[i][0]
      attr[i] = zipped[i][1]
    }
    return cells
  } else {
    cells.sort(compareCells)
    return cells
  }
}
exports.normalize = normalize

//Removes all duplicate cells in the complex
function unique(cells) {
  if(cells.length === 0) {
    return []
  }
  var ptr = 1
    , len = cells.length
  for(var i=1; i<len; ++i) {
    var a = cells[i]
    if(compareCells(a, cells[i-1])) {
      if(i === ptr) {
        ptr++
        continue
      }
      cells[ptr++] = a
    }
  }
  cells.length = ptr
  return cells
}
exports.unique = unique;

//Finds a cell in a normalized cell complex
function findCell(cells, c) {
  var lo = 0
    , hi = cells.length-1
    , r  = -1
  while (lo <= hi) {
    var mid = (lo + hi) >> 1
      , s   = compareCells(cells[mid], c)
    if(s <= 0) {
      if(s === 0) {
        r = mid
      }
      lo = mid + 1
    } else if(s > 0) {
      hi = mid - 1
    }
  }
  return r
}
exports.findCell = findCell;

//Builds an index for an n-cell.  This is more general than dual, but less efficient
function incidence(from_cells, to_cells) {
  var index = new Array(from_cells.length)
  for(var i=0, il=index.length; i<il; ++i) {
    index[i] = []
  }
  var b = []
  for(var i=0, n=to_cells.length; i<n; ++i) {
    var c = to_cells[i]
    var cl = c.length
    for(var k=1, kn=(1<<cl); k<kn; ++k) {
      b.length = bits.popCount(k)
      var l = 0
      for(var j=0; j<cl; ++j) {
        if(k & (1<<j)) {
          b[l++] = c[j]
        }
      }
      var idx=findCell(from_cells, b)
      if(idx < 0) {
        continue
      }
      while(true) {
        index[idx++].push(i)
        if(idx >= from_cells.length || compareCells(from_cells[idx], b) !== 0) {
          break
        }
      }
    }
  }
  return index
}
exports.incidence = incidence

//Computes the dual of the mesh.  This is basically an optimized version of buildIndex for the situation where from_cells is just the list of vertices
function dual(cells, vertex_count) {
  if(!vertex_count) {
    return incidence(unique(skeleton(cells, 0)), cells, 0)
  }
  var res = new Array(vertex_count)
  for(var i=0; i<vertex_count; ++i) {
    res[i] = []
  }
  for(var i=0, len=cells.length; i<len; ++i) {
    var c = cells[i]
    for(var j=0, cl=c.length; j<cl; ++j) {
      res[c[j]].push(i)
    }
  }
  return res
}
exports.dual = dual

//Enumerates all cells in the complex
function explode(cells) {
  var result = []
  for(var i=0, il=cells.length; i<il; ++i) {
    var c = cells[i]
      , cl = c.length|0
    for(var j=1, jl=(1<<cl); j<jl; ++j) {
      var b = []
      for(var k=0; k<cl; ++k) {
        if((j >>> k) & 1) {
          b.push(c[k])
        }
      }
      result.push(b)
    }
  }
  return normalize(result)
}
exports.explode = explode

//Enumerates all of the n-cells of a cell complex
function skeleton(cells, n) {
  if(n < 0) {
    return []
  }
  var result = []
    , k0     = (1<<(n+1))-1
  for(var i=0; i<cells.length; ++i) {
    var c = cells[i]
    for(var k=k0; k<(1<<c.length); k=bits.nextCombination(k)) {
      var b = new Array(n+1)
        , l = 0
      for(var j=0; j<c.length; ++j) {
        if(k & (1<<j)) {
          b[l++] = c[j]
        }
      }
      result.push(b)
    }
  }
  return normalize(result)
}
exports.skeleton = skeleton;

//Computes the boundary of all cells, does not remove duplicates
function boundary(cells) {
  var res = []
  for(var i=0,il=cells.length; i<il; ++i) {
    var c = cells[i]
    for(var j=0,cl=c.length; j<cl; ++j) {
      var b = new Array(c.length-1)
      for(var k=0, l=0; k<cl; ++k) {
        if(k !== j) {
          b[l++] = c[k]
        }
      }
      res.push(b)
    }
  }
  return normalize(res)
}
exports.boundary = boundary;

//Computes connected components for a dense cell complex
function connectedComponents_dense(cells, vertex_count) {
  var labels = new UnionFind(vertex_count)
  for(var i=0; i<cells.length; ++i) {
    var c = cells[i]
    for(var j=0; j<c.length; ++j) {
      for(var k=j+1; k<c.length; ++k) {
        labels.link(c[j], c[k])
      }
    }
  }
  var components = []
    , component_labels = labels.ranks
  for(var i=0; i<component_labels.length; ++i) {
    component_labels[i] = -1
  }
  for(var i=0; i<cells.length; ++i) {
    var l = labels.find(cells[i][0])
    if(component_labels[l] < 0) {
      component_labels[l] = components.length
      components.push([cells[i].slice(0)])
    } else {
      components[component_labels[l]].push(cells[i].slice(0))
    }
  }
  return components
}

//Computes connected components for a sparse graph
function connectedComponents_sparse(cells) {
  var vertices  = unique(normalize(skeleton(cells, 0)))
    , labels    = new UnionFind(vertices.length)
  for(var i=0; i<cells.length; ++i) {
    var c = cells[i]
    for(var j=0; j<c.length; ++j) {
      var vj = findCell(vertices, [c[j]])
      for(var k=j+1; k<c.length; ++k) {
        labels.link(vj, findCell(vertices, [c[k]]))
      }
    }
  }
  var components        = []
    , component_labels  = labels.ranks
  for(var i=0; i<component_labels.length; ++i) {
    component_labels[i] = -1
  }
  for(var i=0; i<cells.length; ++i) {
    var l = labels.find(findCell(vertices, [cells[i][0]]));
    if(component_labels[l] < 0) {
      component_labels[l] = components.length
      components.push([cells[i].slice(0)])
    } else {
      components[component_labels[l]].push(cells[i].slice(0))
    }
  }
  return components
}

//Computes connected components for a cell complex
function connectedComponents(cells, vertex_count) {
  if(vertex_count) {
    return connectedComponents_dense(cells, vertex_count)
  }
  return connectedComponents_sparse(cells)
}
exports.connectedComponents = connectedComponents

},{"bit-twiddle":"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/simplicial-complex/node_modules/bit-twiddle/twiddle.js","union-find":"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/simplicial-complex/node_modules/union-find/index.js"}],"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/node_modules/uniq/uniq.js":[function(require,module,exports){
"use strict"

function unique_pred(list, compare) {
  var ptr = 1
    , len = list.length
    , a=list[0], b=list[0]
  for(var i=1; i<len; ++i) {
    b = a
    a = list[i]
    if(compare(a, b)) {
      if(i === ptr) {
        ptr++
        continue
      }
      list[ptr++] = a
    }
  }
  list.length = ptr
  return list
}

function unique_eq(list) {
  var ptr = 1
    , len = list.length
    , a=list[0], b = list[0]
  for(var i=1; i<len; ++i, b=a) {
    b = a
    a = list[i]
    if(a !== b) {
      if(i === ptr) {
        ptr++
        continue
      }
      list[ptr++] = a
    }
  }
  list.length = ptr
  return list
}

function unique(list, compare, sorted) {
  if(list.length === 0) {
    return list
  }
  if(compare) {
    if(!sorted) {
      list.sort(compare)
    }
    return unique_pred(list, compare)
  }
  if(!sorted) {
    list.sort()
  }
  return unique_eq(list)
}

module.exports = unique

},{}],"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/triangulate.js":[function(require,module,exports){
"use strict"

var ch = require("incremental-convex-hull")
var uniq = require("uniq")

module.exports = triangulate

function LiftedPoint(p, i) {
  this.point = p
  this.index = i
}

function compareLifted(a, b) {
  var ap = a.point
  var bp = b.point
  var d = ap.length
  for(var i=0; i<d; ++i) {
    var s = bp[i] - ap[i]
    if(s) {
      return s
    }
  }
  return 0
}

function triangulate1D(n, points, includePointAtInfinity) {
  if(n === 1) {
    if(includePointAtInfinity) {
      return [ [-1, 0] ]
    } else {
      return []
    }
  }
  var lifted = points.map(function(p, i) {
    return [ p[0], i ]
  })
  lifted.sort(function(a,b) {
    return a[0] - b[0]
  })
  var cells = new Array(n - 1)
  for(var i=1; i<n; ++i) {
    var a = lifted[i-1]
    var b = lifted[i]
    cells[i-1] = [ a[1], b[1] ]
  }
  if(includePointAtInfinity) {
    cells.push(
      [ -1, cells[0][1], ],
      [ cells[n-1][1], -1 ])
  }
  return cells
}

function triangulate(points, includePointAtInfinity) {
  var n = points.length
  if(n === 0) {
    return []
  }
  
  var d = points[0].length
  if(d < 1) {
    return []
  }

  //Special case:  For 1D we can just sort the points
  if(d === 1) {
    return triangulate1D(n, points, includePointAtInfinity)
  }
  
  //Lift points, sort
  var lifted = new Array(n)
  var upper = 1.0
  for(var i=0; i<n; ++i) {
    var p = points[i]
    var x = new Array(d+1)
    var l = 0.0
    for(var j=0; j<d; ++j) {
      var v = p[j]
      x[j] = v
      l += v * v
    }
    x[d] = l
    lifted[i] = new LiftedPoint(x, i)
    upper = Math.max(l, upper)
  }
  uniq(lifted, compareLifted)
  
  //Double points
  n = lifted.length

  //Create new list of points
  var dpoints = new Array(n + d + 1)
  var dindex = new Array(n + d + 1)

  //Add steiner points at top
  var u = (d+1) * (d+1) * upper
  var y = new Array(d+1)
  for(var i=0; i<=d; ++i) {
    y[i] = 0.0
  }
  y[d] = u

  dpoints[0] = y.slice()
  dindex[0] = -1

  for(var i=0; i<=d; ++i) {
    var x = y.slice()
    x[i] = 1
    dpoints[i+1] = x
    dindex[i+1] = -1
  }

  //Copy rest of the points over
  for(var i=0; i<n; ++i) {
    var h = lifted[i]
    dpoints[i + d + 1] = h.point
    dindex[i + d + 1] =  h.index
  }

  //Construct convex hull
  var hull = ch(dpoints, false)
  if(includePointAtInfinity) {
    hull = hull.filter(function(cell) {
      var count = 0
      for(var j=0; j<=d; ++j) {
        var v = dindex[cell[j]]
        if(v < 0) {
          if(++count >= 2) {
            return false
          }
        }
        cell[j] = v
      }
      return true
    })
  } else {
    hull = hull.filter(function(cell) {
      for(var i=0; i<=d; ++i) {
        var v = dindex[cell[i]]
        if(v < 0) {
          return false
        }
        cell[i] = v
      }
      return true
    })
  }

  if(d & 1) {
    for(var i=0; i<hull.length; ++i) {
      var h = hull[i]
      var x = h[0]
      h[0] = h[1]
      h[1] = x
    }
  }

  return hull
}
},{"incremental-convex-hull":"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/ich.js","uniq":"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/node_modules/uniq/uniq.js"}],"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/parse-svg-path/index.js":[function(require,module,exports){

module.exports = parse

/**
 * expected argument lengths
 * @type {Object}
 */

var length = {a: 7, c: 6, h: 1, l: 2, m: 2, q: 4, s: 4, t: 2, v: 1, z: 0}

/**
 * segment pattern
 * @type {RegExp}
 */

var segment = /([astvzqmhlc])([^astvzqmhlc]*)/ig

/**
 * parse an svg path data string. Generates an Array
 * of commands where each command is an Array of the
 * form `[command, arg1, arg2, ...]`
 *
 * @param {String} path
 * @return {Array}
 */

function parse(path) {
	var data = []
	path.replace(segment, function(_, command, args){
		var type = command.toLowerCase()
		args = parseValues(args)

		// overloaded moveTo
		if (type == 'm' && args.length > 2) {
			data.push([command].concat(args.splice(0, 2)))
			type = 'l'
			command = command == 'm' ? 'l' : 'L'
		}

		while (true) {
			if (args.length == length[type]) {
				args.unshift(command)
				return data.push(args)
			}
			if (args.length < length[type]) throw new Error('malformed path data')
			data.push([command].concat(args.splice(0, length[type])))
		}
	})
	return data
}

function parseValues(args){
	args = args.match(/-?[.0-9]+(?:e[-+]?\d+)?/ig)
	return args ? args.map(Number) : []
}

},{}],"/Users/Romain/Documents/Programmation/ants/AntColony/src/ant.js":[function(require,module,exports){
'use strict';

var sign = require('./utilities.js').sign;
var calculateDistance = require('./utilities.js').distance;

// var points = require('./initializePoints.js').points;
// var citySet = require('./initializePoints.js').citySet;
// var textPointsId = require('./initializePoints.js').textPointsId;
// var possibleStartPointsId = require('./initializePoints.js').possibleStartPointsId;

var liveMousePosition = require('./mouse.js');

var Vector = require('./vector.js');

var random = Math.random;
var floor = Math.floor;
// var REPULSION = 0.05;
// var REPULSIONSPEED = 0.002;
// var ANTVELOCITY = 0.001;

module.exports = function(container, initPoints, options){

    console.log('Options ant :', options);
    // Define those parameters as attributes of Ant object ?
    var REPULSION = options.repSize;
    var REPULSIONSPEED = options.repSpeed;
    var ANTVELOCITY = options.velocity;
    var WEIGHT = options.weight;

    var mouse = liveMousePosition(container);

    var points = initPoints.points;
    var citySet = initPoints.citySet;
    var textPointsId = initPoints.textPointsId;
    var possibleStartPointsId = initPoints.possibleStartPointsId;


    function Ant(point) {
        this.x = point.x;                
        this.y = point.y;
        this.velocity = ANTVELOCITY;
        this.weight = WEIGHT;
        this.repSize = REPULSION;
        this.repSpeed = REPULSIONSPEED;
        this.edge = undefined;
        this.state = "forage";
        this.edges = [];
        this.lastCity = undefined;
        this.origin = point;
        this.destination = undefined;
        this.orientation = undefined;
        this.direction = new Vector(0,0);
        this.prog = 0;
    }
    // forage: the ant wanders around without any pheromon deposition
    // once it finds a city, it starts remembering the nodes it goes through
    // when it finds another city, it computes the path length and adds pheromons one each edges
    // proportionnaly to the shortestness of the path
    // it resets the list of nodes and continues
    // while foraging the ant choses the path with a pheromon preference


    // static methods
    Ant.generateRandStartPoint = function() {
        var randId = Math.floor(possibleStartPointsId.length * random());
        var randStartPoint = points[possibleStartPointsId[randId]];
        return randStartPoint;
    }


    // methods
    Ant.prototype = {

        transit: function(){
            switch (this.state) {
            case "forage":
                var res = this.move();
                if (res.cityReached) {
                    this.state = "pheromoning";
                    this.lastCity = this.origin.id;
                };
                break;
            case "pheromoning":
                var res = this.move();
                if (res.edgeChanged) {
                    this.edges.push(this.edge);
                    // found a city
                    if (res.cityReached && (this.origin.id != this.lastCity) ){
                        // compute the length of the path
                        var pathLength = this.edges.map(function(e){return e.distance}).reduce(function(a,b){return a + b});
                        var deltaPheromone = 1/pathLength;
                        var antWeight = this.weight;
                        this.edges.forEach(function(e){
                            var a = e.pt1, b = e.pt2, weight = 1;  
                            // increased dropped pheromons for textEdges
                            if ((citySet.indexOf(a.id) != -1) && citySet.indexOf(b.id) != -1 && (Math.abs(a.id - b.id) == 1))
                            {
                                weight *= antWeight;
                            }
                            e.pheromon += (deltaPheromone * weight);
                        });

                        this.edges = [this.edge];
                        this.lastCity = this.origin.id;
                    }
                }
              break;
            }

        },

        setDirection: function(){
            var possibleEdges = [];

            for (var i = 0; i < this.origin.nexts.length; i++)
            {
                possibleEdges[i] = this.origin.nexts[i];
            } 

            // console.log('smells1: ', possibleEdges);

            possibleEdges.splice(possibleEdges.indexOf(this.edge),1);

            // flip a coin and either take the smelliest path or a random one
            if (random() < 0.8){
                var smells = possibleEdges.map(function(e){return e.pheromon;});
                var index = smells.indexOf(Math.max.apply(Math, smells));
                this.edge = possibleEdges[index];
            } 
            else{
                this.edge = possibleEdges[floor(random()*possibleEdges.length)];
            }
                

            // set the destination point, being edge.pt1 or edge.pt2
            this.destination = (this.origin == this.edge.pt1) ? this.edge.pt2 : this.edge.pt1;

            this.direction.x = this.destination.x - this.origin.x; 
            this.direction.y = this.destination.y - this.origin.y;

            this.direction.normalize();
        },

        move: function(){
            // console.log('move');
            var edgeChanged;
            var cityReached = false;

            this.direction.x = this.destination.x - this.x; 
            this.direction.y = this.destination.y - this.y;
            this.direction.normalize();

            // on edge
            if ((calculateDistance(this, this.destination) > this.repSpeed)){

                // a delta movement will be applied if collision with obstacle detected
                var delta = this.avoidObstacle();

                this.x += this.velocity * this.direction.x + delta.x * this.repSpeed;
                this.y += this.velocity * this.direction.y + delta.y * this.repSpeed;

                this.prog = this.calculateProgression();
                
                edgeChanged = false;

            // on vertex
            } else {
                // console.log('reached');
                this.step = 0;
                this.prog = 0;
                this.origin = this.destination;
                this.x = this.origin.x;
                this.y = this.origin.y;

                this.setDirection();

                cityReached = (citySet.indexOf(this.origin.id) != -1);
                edgeChanged = true;
            }
            return {cityReached: cityReached, edgeChanged: edgeChanged};
        },

        avoidObstacle: function(){
            var distance = calculateDistance(this, mouse);
        
            if (distance <= this.repSize) {

                return {
                    // delta movement is composed of a repulsion delta and a circular delta 
                    x: (this.x - mouse.x)/distance + (this.y - mouse.y)/distance * 1,
                    y: (this.y - mouse.y)/distance - (this.x - mouse.x)/distance * 1
                };
            }
            else
                return {x:0, y:0};
        },

        calculateProgression: function(){
            var v = new Vector(this.x - this.origin.x, this.y - this.origin.y);
            var norm = v.norm();

            var theta = (v.x * this.edge.direction.x + v.y * this.edge.direction.y) / norm;
            var prog = norm * Math.abs(theta);
            // returns length of projection on edge
            return prog;
        }

    };
    return Ant;
}


},{"./mouse.js":"/Users/Romain/Documents/Programmation/ants/AntColony/src/mouse.js","./utilities.js":"/Users/Romain/Documents/Programmation/ants/AntColony/src/utilities.js","./vector.js":"/Users/Romain/Documents/Programmation/ants/AntColony/src/vector.js"}],"/Users/Romain/Documents/Programmation/ants/AntColony/src/antsGroup.js":[function(require,module,exports){
'use strict'

module.exports = function (Ant) {

	var nbAntsPerStep = 100;

	function createGroup(population){
		for (var i = 0; i < nbAntsPerStep; i++) {
			var newAnt = new Ant(Ant.generateRandStartPoint());
			newAnt.setDirection();
			population.push(newAnt);
		}

// 		console.log('Created Ants Group: \
// (+ ' + nbAntsPerStep + ') => ' + population.length);

		return population;
	}

	function removeGroup(population, nbDead){
		population = population.slice(0, population.length - nbDead);

// 		console.log('Removed Ants Group: \
// (- ' + nbAntsPerStep + ') => ' + population.length);

		return population;

	}

	return {
		create: createGroup,
		remove: removeGroup
	};

}
	
},{}],"/Users/Romain/Documents/Programmation/ants/AntColony/src/createEdges.js":[function(require,module,exports){
'use strict'

var dt = require("delaunay-triangulate");

var range = require('./utilities.js').range;

// var points = require('./initializePoints.js').points;
var textMesh = require('./initializePoints.js').textMesh;
var citySet = require('./initializePoints.js').citySet;
var nbRandomPoints = require('./initializePoints.js').nbRandomPoints;
var forcedEdges = require('./initializePoints.js').forcedEdges;

var Edge = require('./edge.js');


module.exports = function(points){
    // triangulate
    var cells = dt(points.map(function(p){
        return [p.x, p.y];
    }));

    var edges = [];
    var permutations = [[0,1], [0,2], [1,2]];

    // force the edges of the text to be edges of the graph
    if (textMesh) {
        range(0, points.length - nbRandomPoints).forEach(function(id){
            var directLink = forcedEdges[id];
            var textEdge = Edge.create(points[id], points[directLink]);
            edges.push(textEdge);
            points[id].nexts.push(textEdge);
        })
    }


    cells.forEach(function(cell){
       
        for (var i = 0; i < 3; ++i){  // for each point.id listed in current cell
            var pt = points[cell[i]];

            for (var j = 1; j < 3; ++j){ 

                var ptj = points[cell[( i + j ) % 3]]; // pick one of the other 2 points of the cell
                var newEdge = undefined;

                // if pt already has nextEdges ...
                if (pt.nexts.length != 0) {
                    
                    // ... get the points corresponding ...
                    var tempPoints = pt.nexts.map(function(e){
                        return [e.pt1, e.pt2];
                    }).reduce(function(a, b){
                         return a.concat(b);
                    });

                    // ... and check if ptj already is part of the existing nextEdges. If not, add the edge.
                    if (tempPoints.indexOf(ptj) == -1){
                        newEdge = Edge.create(pt, ptj);
                        edges.push(newEdge);
                        pt.nexts.push(newEdge);
                    }
                }
                else {
                    newEdge = Edge.create(pt, ptj);
                    edges.push(newEdge);
                    pt.nexts.push(newEdge);
                }

                // add also the edge to the edge's other point's nextEdges
                if (newEdge != undefined){
                    ptj.nexts.push(newEdge);
                }         
            }

            // add the textEdges to nextEdges map
            if (textMesh && (citySet.indexOf(pt) != -1)) {
                var textEdge = Edge.create(pt, points[pt.id + 1]);
                edges.push(textEdge);
                pt.nexts.push(textEdge);
            }

        }
    });

    return edges;
};
},{"./edge.js":"/Users/Romain/Documents/Programmation/ants/AntColony/src/edge.js","./initializePoints.js":"/Users/Romain/Documents/Programmation/ants/AntColony/src/initializePoints.js","./utilities.js":"/Users/Romain/Documents/Programmation/ants/AntColony/src/utilities.js","delaunay-triangulate":"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/triangulate.js"}],"/Users/Romain/Documents/Programmation/ants/AntColony/src/edge.js":[function(require,module,exports){
'use strict';

var sqrt = Math.sqrt;
var pow = Math.pow;
var abs = Math.abs;
var atan = Math.atan;

var Vector = require('./vector.js');


function Edge(ptA, ptB) {
    var distance = sqrt( pow(ptA.x - ptB.x, 2) + pow(ptA.y - ptB.y, 2) );

    // find line equation ax + by + c = 0
    var a = 1;
    var b = - (ptB.x - ptA.x) / (ptB.y - ptA.y);

    // orientate vector (a,b)
    if (b < 0){
        b = -b;
        a = -a;
    }

    // normalize vector (a,b)
    var n = new Vector(a, b);
    n.normalize();

    var c = - (a * ptA.x + b * ptA.y);

    // // calculate vector director
    var v = new Vector(ptB.x - ptA.x, ptB.y - ptA.y);
    
    v.normalize();

    this.id = undefined;
    this.pt1 = ptA;
    this.pt2 = ptB;
    this.direction = v;
    this.orthDirection = n; 
    this.distance = distance;
    this.pheromon = 1/distance;
    this.line = {
        a: a,
        b: b,
        c: c,
    };

    if (this.distance === 0) console.log('ZERO !');
}


// static methods
Edge.create = function(ptA, ptB) {
    var edge = new Edge(ptA, ptB);
    return edge;
}


// methods
Edge.prototype = {

    getOtherPoint: function(point) {
        if (point == this.pt1)
            return this.pt2;
        else if (point == this.pt2)
            return this.pt1;
        else
            console.log("Error");
    },

    calculateDistance: function(x, y) {
        var a = this.line.a,
            b = this.line.b,
            c = this.line.c;
        return abs(a * x + b * y + c) / Math.sqrt(Math.pow(a,2) + Math.pow(b,2));
    },

}
module.exports = Edge;
},{"./vector.js":"/Users/Romain/Documents/Programmation/ants/AntColony/src/vector.js"}],"/Users/Romain/Documents/Programmation/ants/AntColony/src/initializePoints.js":[function(require,module,exports){
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

},{"./point.js":"/Users/Romain/Documents/Programmation/ants/AntColony/src/point.js","./svgPath.js":"/Users/Romain/Documents/Programmation/ants/AntColony/src/svgPath.js","./utilities.js":"/Users/Romain/Documents/Programmation/ants/AntColony/src/utilities.js","parse-svg-path":"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/parse-svg-path/index.js"}],"/Users/Romain/Documents/Programmation/ants/AntColony/src/mouse.js":[function(require,module,exports){
'use strict'

module.exports = function (container){

	var mouse = {
	    x: 0,
	    y: 0
	};

	container.addEventListener( 'mousemove', function(e){
	    var rect = container.getBoundingClientRect();

	    mouse.x = (e.clientX - rect.left ) / rect.width;
	    mouse.y = (e.clientY - rect.top )/ rect.height;
	});

	return mouse;

};

},{}],"/Users/Romain/Documents/Programmation/ants/AntColony/src/point.js":[function(require,module,exports){
'use strict'

function Point(x, y) {
    this.id = undefined;                
    this.x = x;
    this.y = y;
    this.nexts = [];
}

module.exports = Point;
},{}],"/Users/Romain/Documents/Programmation/ants/AntColony/src/rendering.js":[function(require,module,exports){
'use strict'

var antFunction = require('./ant.js');
var antsGroup = require('./antsGroup');

var random = Math.random;

var RANDOMMVT = 0.003;
var ANTSIZE = 0.002;

module.exports = function(container, pointsMap, options){

	if(!container)
		throw new TypeError('Missing container');

	// Ants variables
	var edges = pointsMap.edges;
	var objPopulationInitial = options.nbAnts;
	var objPopulation = objPopulationInitial;
	var pointsInfos = pointsMap.pointsInfos;
	var population = [];
	var nbAntsPerStep = 100;
	
	var Ant = antFunction(container, pointsInfos, options);
	antsGroup = antsGroup(Ant);

	// Animation variables
	var animID;
	var deltaTime;
	var FPSCount;
	var lastUpdate = performance.now();
	// var FPSMonitor = document.querySelector('#FPS');
	// var dTMonitor = document.querySelector('#dT');
	var refreshTime = 0;
	var maxDeltaTime = 40;
	var FPSOverLimitCount = 0;
	var FPSUnderLimitCount = 0;


	// Canvas
	var canvasList = document.getElementsByTagName("canvas");
	
	if (canvasList.length === 0){
		var canvas = document.createElement("canvas");
		var rect = container.getBoundingClientRect();
		canvas.width = rect.width;
		canvas.height = rect.height;
		canvas.style.backgroundColor = "rgba(250, 250, 250, 0)"; 
		container.appendChild(canvas);
	}
	else{
		var canvas = canvasList[0];
		console.log('CANVAS');
	}
	var context = canvas.getContext("2d");
	context.clearRect ( 0 , 0 , canvas.width, canvas.height );
	

	function checkAntNumber(antNumber){
		if (antNumber < objPopulation - 50){
			// FPSMonitor.style.color = "green";
			population = antsGroup.create(population);
		}	
		else if (antNumber > objPopulation){
			population = antsGroup.remove(population, antNumber - objPopulation);
			// FPSMonitor.style.color = "red";
		}
		// else
		// 	FPSMonitor.style.color = "white";
	}

	function displayFPS(dT){
		FPSCount = (1000/dT).toFixed(2);
		var t = dT.toFixed(2);
		// FPSMonitor.textContent = 'FPS : ' + FPSCount;  
		// dTMonitor.textContent = 'nbAnts : ' + population.length;
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
			objPopulation = objPopulation * maxDeltaTime / deltaTime;
			FPSOverLimitCount = 0;
		}

		while (FPSUnderLimitCount > 50 && objPopulation < objPopulationInitial) {
			objPopulation += 10;
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
		context.fillStyle = "rgba(250, 250, 250, 0.4)";
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
		objPopulation = opts.nbAnts;

		population.forEach(function(ant){
			ant.velocity = opts.velocity;
			ant.weight = opts.weight;
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

},{"./ant.js":"/Users/Romain/Documents/Programmation/ants/AntColony/src/ant.js","./antsGroup":"/Users/Romain/Documents/Programmation/ants/AntColony/src/antsGroup.js"}],"/Users/Romain/Documents/Programmation/ants/AntColony/src/svgPath.js":[function(require,module,exports){
'use strict';

var point1 = "m 18.25,19.5 c 18.25,15.7677688 15.0000951,12.75 11,12.75 c 6.99990488,12.75 3.75,15.7677688 3.75,19.5 c 3.75,23.2322312 6.99990488,26.25 11,26.25 c 15.0000951,26.25 18.25,23.2322312 18.25,19.5 z m 4.25,19.5 c 4.25,16.0525294 7.26810863,13.25 11,13.25 c 14.7318914,13.25 17.75,16.0525294 17.75,19.5 c 17.75,22.9474706 14.7318914,25.75 11,25.75 c 7.26810863,25.75 4.25,22.9474706 4.25,19.5 z";
var point2 = "m 89.25,8.5 c 89.25,4.21605787 85.5528716,0.75 81,0.75 c 76.4471284,0.75 72.75,4.21605787 72.75,8.5 c 72.75,12.7839421 76.4471284,16.25 81,16.25 c 85.5528716,16.25 89.25,12.7839421 89.25,8.5 z m 73.25,8.5 c 73.25,4.49967088 76.7163156,1.25 81,1.25 c 85.2836844,1.25 88.75,4.49967088 88.75,8.5 c 88.75,12.5003291 85.2836844,15.75 81,15.75 c 76.7163156,15.75 73.25,12.5003291 73.25,8.5 z";
var point3 = "m 160.25,11 c 160.25,5.33348375 155.208168,0.75 149,0.75 c 142.791832,0.75 137.75,5.33348375 137.75,11 c 137.75,16.6665162 142.791832,21.25 149,21.25 c 155.208168,21.25 160.25,16.6665162 160.25,11 z m 138.25,11 c 138.25,5.62082125 143.057903,1.25 149,1.25 c 154.942097,1.25 159.75,5.62082125 159.75,11 c 159.75,16.3791787 154.942097,20.75 149,20.75 c 143.057903,20.75 138.25,16.3791787 138.25,11 z";
var point4 = "m 160.25,76.5 c 160.25,72.7677688 157.000095,69.75 153,69.75 c 148.999905,69.75 145.75,72.7677688 145.75,76.5 c 145.75,80.2322312 148.999905,83.25 153,83.25 c 157.000095,83.25 160.25,80.2322312 160.25,76.5 z m 146.25,76.5 c 146.25,73.0525294 149.268109,70.25 153,70.25 c 156.731891,70.25 159.75,73.0525294 159.75,76.5 c 159.75,79.9474706 156.731891,82.75 153,82.75 c 149.268109,82.75 146.25,79.9474706 146.25,76.5 z";
var point5 = "m 95.25,76 c 95.25,70.3362795 90.4344065,65.75 84.5,65.75 c 78.5655935,65.75 73.75,70.3362795 73.75,76 c 73.75,81.6637205 78.5655935,86.25 84.5,86.25 c 90.4344065,86.25 95.25,81.6637205 95.25,76 z m 74.25,76 c 74.25,70.6180255 78.8364268,66.25 84.5,66.25 c 90.1635732,66.25 94.75,70.6180255 94.75,76 c 94.75,81.3819745 90.1635732,85.75 84.5,85.75 c 78.8364268,85.75 74.25,81.3819745 74.25,76 z";
var point6 = "m 20.25,75 c 20.25,70.9919338 16.7764995,67.75 12.5,67.75 c 8.22350046,67.75 4.75,70.9919338 4.75,75 c 4.75,79.0080662 8.22350046,82.25 12.5,82.25 c 16.7764995,82.25 20.25,79.0080662 20.25,75 z m 5.25,75 c 5.25,71.2760797 8.49222829,68.25 12.5,68.25 c 16.5077717,68.25 19.75,71.2760797 19.75,75 c 19.75,78.7239203 16.5077717,81.75 12.5,81.75 c 8.49222829,81.75 5.25,78.7239203 5.25,75 z";
var point7 = "m 23.25,139 c 23.25,132.786797 18.2132034,127.75 12,127.75 c 5.78679656,127.75 0.75,132.786797 0.75,139 c 0.75,145.213203 5.78679656,150.25 12,150.25 c 18.2132034,150.25 23.25,145.213203 23.25,139 z m 1.25,139 c 1.25,133.062939 6.06293894,128.25 12,128.25 c 17.9370611,128.25 22.75,133.062939 22.75,139 c 22.75,144.937061 17.9370611,149.75 12,149.75 c 6.06293894,149.75 1.25,144.937061 1.25,139 z";
var point8 = "m 95.2429209,133.141624 c 95.5537927,127.209836 90.7562698,122.141736 84.5337404,121.815627 c 78.311211,121.489518 73.0102086,126.028377 72.6993368,131.960165 c 72.388465,137.891953 77.1859879,142.960053 83.4085173,143.286162 c 89.6310467,143.612271 94.932049,139.073412 95.2429209,133.141624 z m 73.1986516,131.986333 c 73.4947711,126.336036 78.555388,122.003001 84.5075724,122.314942 c 90.4597567,122.626883 95.0397256,127.465159 94.7436061,133.115456 c 94.4474866,138.765754 89.3868696,143.098788 83.4346853,142.786847 c 77.4825009,142.474907 72.902532,137.63663 73.1986516,131.986333 z";
var point9 = "m 167.728023,135.062479 c 168.038327,129.141533 163.482117,124.089864 157.551662,123.779062 c 151.621208,123.468261 146.561914,128.016002 146.25161,133.936948 c 145.941307,139.857894 150.497517,144.909562 156.427971,145.220364 c 162.358426,145.531166 167.41772,140.983425 167.728023,135.062479 z m 146.750925,133.963116 c 147.046767,128.318121 151.870617,123.982018 157.525494,124.278377 c 163.180372,124.574737 167.52455,129.391317 167.228709,135.036311 c 166.932867,140.681305 162.109017,145.017409 156.454139,144.72105 c 150.799262,144.42469 146.455083,139.60811 146.750925,133.963116 z";
var lettresANT = "m 25.6185414,31.9447266 l 25.6185414,31.1003906 l 25.388268,30.4095703 l 24.9277211,30.0257812 l 24.5439321,29.5652344 l 24.0833852,29.1814453 l 23.6995961,28.7208984 l 23.0087758,28.490625 l 22.1453678,28.490625 l 21.1512711,28.490625 l 20.475768,28.490625 l 19.8617055,28.7208984 l 19.4779164,29.1814453 l 18.7870961,29.3349609 l 17.5478579,29.3349609 l 16.669907,29.3349609 l 16.0024971,29.3349609 l 15.3349802,29.3349609 l 14.5654164,29.3349609 l 13.9513539,29.5652344 l 13.5675649,30.0257812 l 12.8767446,30.1792969 l 12.0324086,30.1792969 l 11.4183461,30.4095703 l 11.0345571,30.8701172 l 10.3437367,31.0236328 l 9.72967424,31.2539062 l 9.34588518,31.7144531 l 8.8853383,32.0982422 l 8.50154924,32.5587891 l 8.04100237,32.9425781 l 7.88748674,33.6333984 l 7.6572133,34.2474609 l 7.19666643,34.63125 l 7.0431508,35.3220703 l 7.0431508,36.1664062 l 7.19666643,36.7804687 l 7.6572133,37.1642578 l 8.04100237,37.6248047 l 8.50154924,38.0085937 l 8.8853383,38.4691406 l 9.34588518,38.8529297 l 9.72967424,39.3134766 l 10.3437367,39.4669922 l 11.4725961,39.4669922 l 12.3122961,39.4669922 l 12.8743818,39.4669922 l 13.6496873,39.4669922 l 14.5398528,39.4669922 l 15.4097524,39.4669922 l 16.1005727,39.3134766 l 16.4843617,38.8529297 l 17.0984242,38.6226562 l 17.7892446,38.4691406 l 18.1730336,38.0085937 l 18.6335805,37.6248047 l 19.0173696,37.1642578 l 19.6314321,36.9339844 l 20.3222524,36.7804687 l 20.7060414,36.3199219 l 21.3201039,36.0896484 l 22.0109242,35.9361328 l 22.3947133,35.4755859 l 23.0087758,35.2453125 l 23.6995961,35.0917969 l 23.9298696,34.4777344 l 24.0833852,33.7869141 l 24.5439321,33.403125 l 24.9277211,32.9425781 l 25.388268,32.5587891 l 25.6185414,31.9447266 z m 37.5927602,43.9189453 l 37.2089711,44.3794922 l 36.7484242,44.7632813 l 36.3646352,45.2238281 l 35.6738149,45.3773438 l 35.0597524,45.6076172 l 34.6759633,46.0681641 l 33.985143,46.2216797 l 33.395172,46.2216797 l 32.4874219,46.2216797 l 31.4521352,46.2216797 l 30.8380727,46.0681641 l 30.4542836,45.6076172 l 29.7634633,45.3773438 l 29.1494008,45.2238281 l 28.7656117,44.7632813 l 28.3050649,44.3794922 l 27.9212758,43.9189453 l 27.2304555,43.6886719 l 26.616393,43.5351563 l 26.2326039,43.0746094 l 25.7720571,43.0746094 l 25.388268,43.5351563 l 24.6974477,43.6886719 l 23.8531117,43.6886719 l 23.2390492,43.9189453 l 22.8552602,44.3794922 l 22.1644399,44.5330078 l 21.5503774,44.7632813 l 21.1665883,45.2238281 l 20.475768,45.3773438 l 19.6314321,45.3773438 l 19.0173696,45.6076172 l 18.6335805,46.0681641 l 17.9427602,46.2216797 l 17.1279087,46.2216797 l 16.2540883,46.2216797 l 15.6400258,46.4519531 l 15.2562367,46.9125 l 14.7956899,46.9125 l 14.4119008,46.4519531 l 13.7210805,46.2216797 l 12.7926837,46.2216797 l 12.01698,46.2216797 l 11.36469,46.2216797 l 10.48491,46.2216797 l 9.4994008,46.2216797 l 8.8853383,46.0681641 l 8.50154924,45.6076172 l 7.81072893,45.3773438 l 7.19666643,45.2238281 l 6.81287737,44.7632813 l 6.12205705,44.5330078 l 5.50799455,44.3794922 l 5.12420549,43.9189453 l 4.66365862,43.5351563 l 4.27986955,43.0746094 l 3.81932268,42.6908203 l 3.43553362,42.2302734 l 2.97498674,41.8464844 l 2.59119768,41.3859375 l 2.1306508,41.0021484 l 1.97713518,40.3880859 l 1.74686174,39.6972656 l 1.28631487,39.3134766 l 1.13279924,38.6994141 l 1.13279924,37.8550781 l 0.902525804,37.1642578 l 0.441978929,36.7804688 l 0.288463304,36.1664063 l 0.288463304,35.2521433 l 0.288463304,34.3970875 l 0.288463304,33.6333984 c 0.288463304,33.5008705 0.441978929,32.9425781 0.441978929,32.9425781 l 0.902525804,32.5587891 l 1.13279924,31.9447266 l 1.13279924,31.1003906 l 1.28631487,30.4095703 l 1.74686174,30.0257813 l 1.97713518,29.4117188 l 2.1306508,28.7208984 l 2.59119768,28.3371094 l 2.97498674,27.8765625 l 3.43553362,27.4927734 l 3.81932268,27.0322266 l 4.27986955,26.6484375 l 4.66365862,26.1878906 l 5.27772112,25.9576172 l 6.12205705,25.9576172 l 6.81287737,25.8041016 l 7.19666643,25.3435547 l 7.81072893,25.1132813 l 8.65506487,25.1132813 l 9.34588518,24.9597656 l 9.72967424,24.4992188 l 10.3437367,24.2689453 l 10.9370585,24.2689453 l 11.6695294,24.2689453 l 12.1244404,24.2689453 l 12.7890925,24.2689453 l 13.7210805,24.2689453 c 13.9566084,24.2689453 14.4119008,24.1154297 14.4119008,24.1154297 l 14.7956899,23.6548828 l 15.4097524,23.4246094 l 15.8914929,23.4246094 l 16.3780342,23.4246094 l 16.7165068,23.4246094 l 17.3264618,23.4246094 l 17.9427602,23.4246094 c 18.1782881,23.4246094 18.6335805,23.2710938 18.6335805,23.2710938 l 19.0173696,22.8105469 l 19.6314321,22.5802734 l 20.1258875,22.5802734 l 20.6556616,22.5802734 l 21.3026338,22.5802734 l 22.4916732,22.5802734 l 23.4347831,22.5802734 l 23.8531117,22.5802734 l 24.5439321,22.4267578 l 24.9277211,21.9662109 l 25.388268,21.5824219 l 25.7720571,21.121875 l 26.3861196,20.8916016 l 27.0769399,20.7380859 l 27.3072133,20.1240234 l 27.3072133,19.3201067 l 27.3072133,18.4353516 l 27.0769399,17.7445313 l 26.616393,17.3607422 l 26.616393,16.9001953 l 27.0769399,16.5164063 l 27.0769399,16.0558594 l 26.616393,15.6720703 l 26.4628774,15.0580078 l 26.2326039,14.3671875 l 25.7720571,13.9833984 l 25.388268,13.5228516 l 24.9277211,13.1390625 l 24.5439321,12.6785156 l 23.8531117,12.4482422 l 23.0087758,12.4482422 l 22.3947133,12.2947266 l 22.0109242,11.8341797 l 21.5503774,11.8341797 l 21.1665883,12.2947266 l 20.475768,12.4482422 l 19.5356151,12.4482422 l 18.6789018,12.4482422 l 18.1006461,12.4482422 l 17.4906911,12.4482422 l 16.8692027,12.4482422 l 16.2962649,12.4482422 l 15.7558552,12.4482422 l 15.1644088,12.4482422 l 14.5654164,12.4482422 l 13.9513539,12.6785156 l 13.5675649,13.1390625 l 13.107018,13.5228516 l 12.7232289,13.9833984 l 12.2626821,14.3671875 l 12.1091664,15.0580078 l 12.1091664,15.7020899 l 12.1091664,16.694579 l 12.1091664,17.5910156 l 11.878893,18.2050781 l 11.4183461,18.5888672 l 11.0345571,19.0494141 l 10.5740102,19.4332031 l 10.1902211,19.89375 l 9.4994008,20.0472656 l 8.6444784,20.0472656 l 7.81072893,20.0472656 l 7.19666643,19.89375 l 6.81287737,19.4332031 l 6.12205705,19.2029297 l 5.50799455,19.0494141 l 5.12420549,18.5888672 l 4.66365862,18.2050781 l 4.51014299,17.5910156 l 4.51014299,16.7466797 l 4.27986955,16.0558594 l 3.81932268,15.6720703 l 3.66580705,15.0580078 l 3.66580705,14.2136719 l 3.81932268,13.5228516 l 4.27986955,13.1390625 l 4.51014299,12.525 l 4.66365862,11.8341797 l 5.12420549,11.4503906 l 5.50799455,10.9898438 l 5.96854143,10.6060547 l 6.35233049,10.1455078 l 6.96639299,9.91523438 l 7.6572133,9.76171875 l 8.04100237,9.30117188 l 8.65506487,9.07089844 l 9.34588518,8.91738281 l 9.72967424,8.45683594 l 10.3437367,8.2265625 l 10.9890621,8.2265625 l 11.4361001,8.2265625 l 11.9428075,8.2265625 l 12.2681583,8.2265625 l 12.9549793,8.2265625 l 13.7210805,8.2265625 l 14.4119008,8.07304688 l 14.7956899,7.6125 l 15.2562367,7.6125 l 15.6400258,8.07304688 l 16.2540883,8.2265625 l 16.960917,8.2265625 l 17.4915199,8.2265625 l 18.029098,8.2265625 l 18.7932683,8.2265625 l 19.6314321,8.2265625 l 20.3222524,8.45683594 l 20.7060414,8.91738281 l 21.3201039,9.07089844 l 22.0440828,9.07089844 l 22.879456,9.07089844 l 23.8531117,9.07089844 c 24.0950474,9.07089844 24.5439321,9.30117188 24.5439321,9.30117188 l 24.9277211,9.76171875 l 25.5417836,9.91523438 l 26.2326039,10.1455078 l 26.616393,10.6060547 l 27.0769399,10.9898438 l 27.4607289,11.4503906 l 27.9212758,11.8341797 l 28.3050649,12.2947266 l 28.7656117,12.6785156 l 29.1494008,13.1390625 l 29.6099477,13.5228516 l 29.8402211,14.2136719 l 29.9937367,14.8277344 l 30.4542836,15.2115234 l 30.6845571,15.9023438 l 30.6845571,16.5032774 c 30.6845571,16.7271532 30.6845571,16.951029 30.6845571,17.1749048 l 30.6845571,17.8765741 l 30.6845571,18.5420549 l 30.6845571,19.1652008 c 30.6845571,19.7091147 30.6845571,20.1240234 30.6845571,20.1240234 l 30.8380727,20.7380859 l 31.2986196,21.121875 l 31.528893,21.8126953 l 31.528893,22.6060737 c 31.528893,22.862869 31.528893,23.1196643 31.528893,23.3764596 l 31.528893,24.110588 c 31.528893,24.3770742 31.528893,24.6555846 31.528893,24.943475 l 31.528893,25.8874136 l 31.528893,26.7741689 c 31.528893,27.0744041 31.528893,27.3746394 31.528893,27.6748747 l 31.528893,28.5871138 l 31.528893,29.2621252 l 31.528893,30.0192512 c 31.528893,30.2644768 31.528893,30.5128063 31.528893,30.7612527 l 31.528893,31.5916534 l 31.528893,32.2756429 l 31.528893,33.1500361 l 31.528893,33.8138594 l 31.528893,34.5533746 c 31.528893,35.0310778 31.528893,35.3220703 31.528893,35.3220703 l 31.6824086,35.9361328 l 32.1429555,36.3199219 l 32.5267446,36.7804688 l 32.9872914,37.1642578 l 33.3710805,37.6248047 l 33.985143,37.7783203 l 34.6759633,38.0085938 l 35.0597524,38.4691406 l 35.6738149,38.6226563 l 36.3646352,38.8529297 l 36.7484242,39.3134766 l 37.3624867,39.4669922 l 38.0533071,39.6972656 l 38.2835805,40.3880859 l 38.2835805,41.0212323 l 38.2835805,42.0753937 l 38.2835805,42.9210938 l 38.0533071,43.5351563 l 37.5927602,43.9189453 z m 87.4181857,40.5416016 l 87.2646701,41.2324219 l 87.0343966,41.8464844 l 86.5738497,42.2302734 l 86.1900607,42.6908203 l 85.7295138,43.0746094 l 85.3457247,43.5351563 l 84.6549044,43.6886719 l 84.0408419,43.9189453 l 83.6570529,44.3794922 l 82.9662326,44.5330078 l 82.3828191,44.5330078 l 81.7615463,44.5330078 l 81.1493993,44.5330078 l 80.3991122,44.5330078 l 79.7698767,44.5330078 l 79.2674188,44.5330078 l 78.5992642,44.5330078 l 77.9471246,44.5330078 l 77.2597342,44.5330078 l 76.6694177,44.5330078 l 76.158997,44.5330078 l 75.584606,44.5330078 l 75.0341926,44.5330078 l 74.507757,44.5330078 l 74.2712009,44.5330078 l 73.7096039,44.5330078 l 72.8342013,44.5330078 c 72.868641,44.5000413 72.2201388,44.3794922 72.2201388,44.3794922 l 71.8363497,43.9189453 l 71.3758029,43.5351563 l 71.2222872,42.9210938 l 70.9920138,42.2302734 l 70.5314669,41.8464844 l 70.3779513,41.2324219 l 70.5314669,40.5416016 l 70.9920138,40.1578125 l 71.3758029,39.6972656 l 71.8363497,39.3134766 l 72.2201388,38.8529297 l 72.6806857,38.4691406 l 73.0644747,38.0085938 l 73.6785372,37.7783203 l 74.3693576,37.6248047 l 74.7531466,37.1642578 l 75.2136935,36.7804688 l 75.5974826,36.3199219 l 76.0580294,35.9361328 l 76.2883029,35.3220703 l 76.2883029,34.4777344 l 76.4418185,33.7869141 l 76.9023654,33.403125 l 77.1326388,32.7890625 l 77.1326388,32.0586517 l 77.1326388,31.1347947 l 77.1326388,30.0680557 l 77.1326388,29.2630136 l 77.1326388,28.3870225 l 77.1326388,27.5049475 l 77.1326388,26.6028315 l 77.1326388,25.7716644 l 77.1326388,24.8237403 l 77.1326388,23.9217138 l 77.1326388,23.0904572 l 77.1326388,22.1313494 l 77.1326388,21.1642789 l 77.1326388,20.2382747 l 77.1326388,19.5429215 l 77.1326388,18.9216488 l 77.1326388,18.1794139 l 77.1326388,17.5910156 l 77.2861544,16.9001953 l 77.7467013,16.5164063 l 77.9769747,15.9023438 l 77.7467013,15.2115234 l 77.2861544,14.8277344 l 77.1326388,14.2136719 l 76.9023654,13.5228516 l 76.4418185,13.1390625 l 76.0580294,12.6785156 l 75.3672091,12.4482422 l 74.7531466,12.2947266 l 74.3693576,11.8341797 l 73.6785372,11.6039063 l 73.0644747,11.4503906 l 72.6806857,10.9898438 l 71.9898654,10.7595703 l 71.2135081,10.7595703 l 70.288488,10.7595703 l 69.7711782,10.7595703 l 69.0539946,10.7595703 l 68.2968184,10.7595703 l 67.5658566,10.7595703 l 66.9045912,10.7595703 l 66.0795138,10.7595703 c 66.1163893,10.810201 65.4654513,10.9898438 65.4654513,10.9898438 l 65.0816622,11.4503906 l 64.3908419,11.6039063 l 63.7767794,11.8341797 l 63.6232638,12.525 l 63.3929904,13.1390625 l 62.9324435,13.5228516 l 62.5486544,13.9833984 l 62.0881076,14.3671875 l 61.7043185,14.8277344 l 61.2437716,15.2115234 l 61.090256,15.9023438 l 60.8599826,16.5164063 l 60.3994357,16.9001953 l 60.0156466,17.3607422 l 59.5550997,17.7445313 l 59.4015841,18.4353516 l 59.1713107,19.0494141 l 58.7107638,19.4332031 l 58.5572482,20.1240234 l 58.3269747,20.7380859 l 57.8664279,21.121875 l 57.7129122,21.8126953 l 57.4826388,22.4267578 l 57.0220919,22.8105469 l 56.8685763,23.5013672 l 56.8685763,23.9897102 l 56.8685763,24.698931 l 56.8685763,25.6548179 l 56.8685763,26.3251197 l 56.8685763,27.0902587 l 56.8685763,27.8004637 l 56.8685763,28.5673828 l 57.0220919,29.1814453 l 57.4826388,29.5652344 l 57.7129122,30.2560547 l 57.7129122,31.1003906 l 57.4826388,31.7144531 l 57.0220919,32.0982422 l 56.8685763,32.7890625 l 57.0220919,33.403125 l 57.4826388,33.7869141 l 57.7129122,34.4777344 l 57.7129122,35.3220703 l 57.8664279,35.9361328 l 58.3269747,36.3199219 l 58.7107638,36.7804688 l 59.1713107,37.1642578 l 59.5550997,37.6248047 l 60.1691622,37.7783203 l 60.8599826,38.0085938 l 61.2437716,38.4691406 l 61.7043185,38.8529297 l 62.0881076,39.3134766 l 62.5486544,39.6972656 l 62.7789279,40.3880859 l 62.9324435,41.0021484 l 63.3929904,41.3859375 l 63.3929904,41.8464844 l 62.9324435,42.2302734 l 62.7789279,42.9210938 l 62.5486544,43.5351563 l 61.8578341,43.6886719 l 61.2437716,43.9189453 l 60.8599826,44.3794922 l 60.1691622,44.5330078 l 59.3248263,44.5330078 l 58.7107638,44.7632813 l 58.3269747,45.2238281 l 57.8664279,45.2238281 l 57.4826388,44.7632813 l 57.0220919,44.7632813 l 56.6383029,45.2238281 l 55.9474826,45.3773438 l 55.1031466,45.3773438 l 54.4890841,45.2238281 l 54.1052951,44.7632813 l 53.6447482,44.7632813 l 53.2609591,45.2238281 l 52.5701388,45.3773438 l 51.4974911,45.3773438 l 50.5906332,45.3773438 l 49.9853754,45.3773438 l 49.1927951,45.3773438 c 49.2268402,45.3443039 48.5787326,45.2238281 48.5787326,45.2238281 l 48.1949435,44.7632813 l 47.5041232,44.5330078 l 46.8900607,44.3794922 l 46.5062716,43.9189453 l 46.0457247,43.5351563 l 45.6619357,43.0746094 l 45.2013888,42.6908203 l 45.0478732,42.0767578 l 45.0478732,41.2324219 l 45.2013888,40.5416016 l 45.6619357,40.1578125 l 46.0457247,39.6972656 l 46.5062716,39.3134766 l 46.8900607,38.8529297 l 47.5041232,38.6226563 l 48.1949435,38.4691406 l 48.5787326,38.0085938 l 49.0392794,37.6248047 l 49.4230685,37.1642578 l 49.8836154,36.7804688 l 50.2674044,36.3199219 l 50.7279513,35.9361328 l 50.9582247,35.3220703 l 50.9582247,34.6998663 l 50.9582247,33.9150441 l 50.9582247,32.9400109 l 50.9582247,32.3075544 l 50.9582247,31.4683352 l 50.9582247,30.5891232 l 50.9582247,29.8138742 l 50.9582247,29.023684 l 50.9582247,28.1238047 l 50.9582247,27.0039211 l 50.9582247,26.034375 c 50.9582249,25.8460125 51.1117404,25.3435547 51.1117404,25.3435547 l 51.5722872,24.9597656 l 51.8025607,24.3457031 l 51.8025607,23.5013672 l 51.5722872,22.8105469 l 51.1117404,22.4267578 l 50.9582247,21.8126953 l 50.9582247,21.0290912 l 50.9582247,20.3358854 l 50.9582247,19.6347167 l 50.9582247,18.9575258 l 50.9582247,18.0005653 l 50.9582247,17.3553148 l 50.9582247,16.7466797 l 50.7279513,16.0558594 l 50.2674044,15.6720703 l 50.1138888,15.0580078 l 49.8836154,14.3671875 l 49.4230685,13.9833984 l 49.0392794,13.5228516 l 48.5787326,13.1390625 l 48.1949435,12.6785156 l 47.7343966,12.2947266 l 47.3506076,11.8341797 l 46.8900607,11.4503906 l 46.5062716,10.9898438 l 45.8154513,10.7595703 l 45.2013888,10.6060547 l 45.2013888,10.1455078 l 45.6619357,9.76171875 l 45.8922091,9.14765625 l 45.8922091,8.30332031 l 46.0457247,7.6125 l 46.5062716,7.22871094 l 46.8900607,6.76816406 l 47.3506076,6.384375 l 47.7343966,5.92382813 l 48.3484591,5.69355469 l 49.0392794,5.54003906 l 49.4230685,5.07949219 l 50.037131,4.84921875 l 50.7279513,5.07949219 l 51.1117404,5.54003906 l 51.7258029,5.69355469 l 52.5701388,5.69355469 l 53.2609591,5.92382813 l 53.6447482,6.384375 l 54.2588107,6.53789063 l 54.949631,6.76816406 l 55.3334201,7.22871094 l 55.7939669,7.6125 l 56.177756,8.07304688 l 56.6383029,8.45683594 l 57.0220919,8.91738281 l 57.4826388,9.30117188 l 57.8664279,9.76171875 l 58.3269747,10.1455078 l 58.7107638,10.6060547 l 59.1713107,10.6060547 l 59.5550997,10.1455078 l 60.1691622,9.91523438 l 61.0134982,9.91523438 l 61.7043185,9.76171875 l 62.0881076,9.30117188 l 62.5486544,8.91738281 l 62.9324435,8.45683594 l 63.3929904,8.07304688 l 63.7767794,7.6125 l 64.3908419,7.38222656 l 65.0816622,7.22871094 l 65.4654513,6.76816406 l 66.0795138,6.53789063 l 66.6456679,6.53789063 l 67.302102,6.53789063 l 68.0591887,6.53789063 l 68.6125216,6.53789063 l 69.3033419,6.384375 l 69.687131,5.92382813 l 70.3011935,5.69355469 l 70.9908197,5.69355469 l 71.5481222,5.69355469 l 72.0436015,5.69355469 l 72.7858364,5.69355469 l 73.6785372,5.69355469 l 74.3693576,5.92382813 l 74.7531466,6.384375 l 75.3672091,6.53789063 l 76.0580294,6.76816406 l 76.4418185,7.22871094 l 77.055881,7.38222656 l 77.7467013,7.6125 l 78.1304904,8.07304688 l 78.5910372,8.45683594 l 78.9748263,8.91738281 l 79.4353732,9.30117188 l 79.8191622,9.76171875 l 80.2797091,10.1455078 l 80.5099826,10.8363281 l 80.6634982,11.4503906 l 81.1240451,11.8341797 l 81.3543185,12.525 l 81.5078341,13.1390625 l 81.968381,13.5228516 l 82.1986544,14.2136719 l 82.1986544,14.9400989 l 82.1986544,15.9791025 l 82.1986544,16.7466797 l 81.968381,17.3607422 l 81.5078341,17.7445313 l 81.5078341,18.2050781 l 81.968381,18.5888672 l 81.968381,19.0494141 l 81.5078341,19.4332031 l 81.3543185,20.1240234 l 81.3543185,20.9653893 l 81.3543185,21.8126953 c 81.3543193,22.4267582 81.5078341,22.4267578 81.5078341,22.4267578 l 81.968381,22.8105469 l 82.1986544,23.5013672 l 82.1986544,24.3423076 l 82.1986544,25.2192828 l 82.1986544,25.9295772 l 82.1986544,27.0713807 l 82.1986544,27.6846908 l 82.1986544,28.5673828 c 82.1986529,29.1814445 82.3521701,29.1814453 82.3521701,29.1814453 l 82.8127169,29.5652344 l 83.0429904,30.2560547 l 83.0429904,31.1003906 l 82.8127169,31.7144531 l 82.3521701,32.0982422 l 82.1986544,32.7890625 l 82.3521701,33.403125 l 82.8127169,33.7869141 l 83.0429904,34.4777344 l 83.196506,35.0917969 l 83.6570529,35.4755859 l 84.0408419,35.9361328 l 84.6549044,36.0896484 l 85.3457247,36.3199219 l 85.7295138,36.7804688 l 86.3435763,36.9339844 l 87.0343966,37.1642578 l 87.2646701,37.8550781 l 87.4181857,38.4691406 l 87.8787326,38.8529297 l 88.109006,39.54375 l 87.8787326,40.1578125 l 87.4181857,40.5416016 z m 130.488924,32.0982422 l 130.335408,32.7890625 l 130.335408,33.6886606 l 130.335408,34.4777344 l 130.105135,35.0917969 l 129.644588,35.4755859 l 129.491072,36.1664063 l 129.491072,36.9742631 l 129.491072,37.8550781 l 129.260799,38.4691406 l 128.800252,38.8529297 l 128.646736,39.54375 l 128.416463,40.1578125 l 127.955916,40.5416016 l 127.572127,41.0021484 l 127.11158,41.3859375 l 126.727791,41.8464844 l 126.267244,42.2302734 l 125.883455,42.6908203 l 125.422908,43.0746094 l 125.039119,43.5351563 l 124.578572,43.9189453 l 124.194783,44.3794922 l 123.503963,44.5330078 l 122.8899,44.7632813 l 122.506111,45.2238281 l 121.815291,45.3773438 l 120.970955,45.3773438 l 120.356892,45.6076172 l 119.973103,46.0681641 l 119.282283,46.2216797 l 118.624185,46.2216797 l 117.941566,46.2216797 l 117.118507,46.2216797 l 116.442043,46.2216797 l 115.707548,46.2216797 l 115.060603,46.2216797 l 114.446541,46.0681641 l 114.062752,45.6076172 l 113.371932,45.3773438 l 112.527596,45.3773438 l 111.913533,45.2238281 l 111.529744,44.7632813 l 110.838924,44.5330078 l 110.224861,44.3794922 l 109.841072,43.9189453 l 109.380525,43.5351563 l 108.996736,43.0746094 l 108.536189,42.6908203 l 108.1524,42.2302734 l 107.691853,41.8464844 l 107.308064,41.3859375 l 106.847517,41.0021484 l 106.463728,40.5416016 l 106.003182,40.1578125 l 105.619392,39.6972656 l 105.158846,39.3134766 l 105.00533,38.6994141 l 104.775057,38.0085938 l 104.31451,37.6248047 l 104.160994,37.0107422 l 104.160994,36.1664063 l 103.930721,35.4755859 l 103.470174,35.0917969 l 103.316658,34.4777344 l 103.316658,33.7440534 l 103.316658,32.7890625 c 103.316658,32.4311707 103.470174,32.0982422 103.470174,32.0982422 l 103.930721,31.7144531 l 104.160994,31.1003906 l 104.160994,30.3962639 l 104.160994,29.5150941 l 104.160994,28.4306574 l 104.160994,27.7230469 c 104.160994,27.5442917 104.31451,27.0322266 104.31451,27.0322266 l 104.775057,26.6484375 l 105.00533,26.034375 l 104.775057,25.3435547 l 104.31451,24.9597656 l 104.160994,24.3457031 l 104.31451,23.6548828 l 104.775057,23.2710938 l 105.00533,22.6570313 l 105.00533,21.5411634 l 105.00533,20.6271416 l 105.00533,19.9922422 l 105.00533,19.2567085 l 105.00533,18.4955965 l 105.00533,17.6477584 l 105.00533,16.8000002 l 105.00533,15.9023438 l 104.775057,15.2115234 l 104.31451,14.8277344 l 104.160994,14.2136719 l 103.930721,13.5228516 l 103.470174,13.1390625 l 103.086385,12.6785156 l 102.395564,12.4482422 l 101.781502,12.2947266 l 101.397713,11.8341797 l 100.706892,11.6039063 l 100.053357,11.6039063 l 99.4288483,11.6039063 l 98.8200067,11.6039063 l 98.1738846,11.6039063 l 97.5598221,11.4503906 l 97.1760331,10.9898438 l 96.7154862,10.6060547 l 96.3316971,10.1455078 l 95.8711503,9.76171875 l 95.4873612,9.30117188 l 95.0268143,8.91738281 l 94.8732987,8.30332031 l 95.0268143,7.6125 l 95.4873612,7.22871094 l 95.8711503,6.76816406 l 96.3316971,6.384375 l 96.7154862,5.92382813 l 97.1760331,5.54003906 l 97.5598221,5.07949219 l 98.1738846,4.84921875 l 98.7337601,4.84921875 l 99.2183875,4.84921875 l 99.7019758,4.84921875 l 100.17853,4.84921875 l 100.604088,4.84921875 l 101.403807,4.84921875 l 102.002657,4.84921875 l 102.395564,4.84921875 l 103.086385,5.07949219 l 103.470174,5.54003906 l 103.930721,5.54003906 l 104.31451,5.07949219 l 104.775057,4.69570313 l 105.00533,4.08164063 l 105.158846,3.39082031 l 105.619392,3.00703125 l 105.849666,2.39296875 l 106.003182,1.70214844 l 106.617244,1.471875 l 107.308064,1.31835938 l 107.538338,0.704296875 l 107.691853,0.0134765625 l 108.1524,0.0134765625 l 108.382674,0.704296875 l 108.536189,1.31835938 l 108.996736,1.70214844 l 109.380525,2.16269531 l 109.841072,2.54648438 l 110.071346,3.23730469 l 110.071346,4.08164063 l 110.224861,4.69570313 l 110.685408,5.07949219 l 111.069197,5.54003906 l 111.68326,5.69355469 l 112.37408,5.92382813 l 112.757869,6.384375 l 113.218416,6.384375 l 113.602205,5.92382813 l 114.216267,5.69355469 l 115.004627,5.69355469 l 115.904939,5.69355469 l 116.59576,5.92382813 l 116.979549,6.384375 l 117.593611,6.53789063 l 118.251782,6.53789063 l 118.732093,6.53789063 l 119.389054,6.53789063 c 119.610092,6.53789063 119.83113,6.53789063 120.052169,6.53789063 c 120.203765,6.53789063 120.381803,6.53789063 120.572446,6.53789063 c 120.764922,6.53789063 120.957398,6.53789063 121.149874,6.53789063 c 121.352135,6.53789063 121.635686,6.53789063 121.902513,6.53789063 l 122.659627,6.53789063 c 122.89757,6.53789063 123.350447,6.76816406 123.350447,6.76816406 l 123.734236,7.22871094 l 124.194783,7.6125 l 124.578572,8.07304688 l 125.039119,8.45683594 l 125.269392,9.14765625 l 125.039119,9.76171875 l 124.578572,10.1455078 l 124.194783,10.6060547 l 123.734236,10.9898438 l 123.350447,11.4503906 l 122.659627,11.6039063 l 121.9882,11.6039063 l 121.349384,11.6039063 l 120.780988,11.6039063 c 120.593601,11.6039063 120.406214,11.6039063 120.218827,11.6039063 l 119.442448,11.6039063 l 118.774137,11.6039063 l 118.176246,11.6039063 l 117.520245,11.6039063 l 116.836667,11.6039063 l 116.238776,11.6039063 c 116.051863,11.6039063 115.868649,11.6039063 115.690922,11.6039063 c 115.490622,11.6039063 115.297292,11.6039063 115.113494,11.6039063 l 114.536066,11.6039063 c 114.374907,11.6039063 114.226441,11.6039063 114.093322,11.6039063 c 113.645967,11.6039063 113.371932,11.6039063 113.371932,11.6039063 l 112.757869,11.8341797 l 112.37408,12.2947266 l 111.913533,12.6785156 l 111.529744,13.1390625 l 111.069197,13.5228516 l 110.685408,13.9833984 l 110.224861,14.3671875 l 110.071346,15.0580078 l 110.071346,15.9902896 l 110.071346,16.7466797 c 110.071345,17.2189228 110.224861,17.3607422 110.224861,17.3607422 l 110.685408,17.7445313 l 110.685408,18.2050781 l 110.224861,18.5888672 l 110.071346,19.2796875 l 110.071346,20.1730483 c 110.071346,20.3950552 110.071346,20.6346766 110.071346,20.8877198 l 110.071346,21.6569848 l 110.071346,22.3905202 l 110.071346,23.1873618 l 110.071346,23.9709346 l 110.071346,24.6535534 l 110.071346,25.4737351 l 110.071346,26.265381 l 110.071346,26.9132294 l 110.071346,27.5442921 c 110.071346,28.1742542 110.071346,28.5673828 110.071346,28.5673828 l 110.224861,29.1814453 l 110.685408,29.5652344 l 110.915682,30.2560547 l 110.915682,31.1003906 l 111.069197,31.7144531 l 111.529744,32.0982422 l 111.760017,32.7890625 l 111.760017,33.6333984 l 111.913533,34.2474609 l 112.37408,34.63125 l 112.604353,35.3220703 l 112.757869,35.9361328 l 113.371932,36.0896484 l 114.062752,36.3199219 l 114.446541,36.7804688 l 115.060603,36.9339844 l 115.742318,36.9339844 l 116.226866,36.9339844 l 116.683917,36.9339844 l 117.041852,36.9339844 l 117.727348,36.9339844 l 118.437947,36.9339844 l 119.128767,36.7804688 l 119.512557,36.3199219 l 119.973103,35.9361328 l 120.356892,35.4755859 l 120.817439,35.0917969 l 121.047713,34.4777344 l 121.201228,33.7869141 l 121.661775,33.403125 l 121.892049,32.7890625 l 121.892049,31.9447266 l 122.045564,31.2539063 l 122.506111,30.8701172 l 122.736385,30.2560547 l 122.8899,29.5652344 l 123.350447,29.1814453 l 123.734236,28.7208984 l 124.348299,28.490625 l 125.039119,28.3371094 l 125.422908,27.8765625 l 126.036971,27.6462891 l 126.881307,27.6462891 l 127.572127,27.8765625 l 127.955916,28.3371094 l 128.569978,28.490625 l 129.260799,28.7208984 l 129.644588,29.1814453 l 130.105135,29.5652344 l 130.335408,30.2560547 l 130.488924,30.8701172 l 130.949471,31.2539063 l 130.949471,31.7144531 l 130.488924,32.0982422 z ";
var lettreS = "m 138.811517,35.1061011 c 138.63781,35.7365111 138.606564,35.8010614 138.198455,36.2091704 c 138.056753,36.3508722 137.999155,37.0215928 137.948798,37.2230199 c 137.869537,37.5400622 137.926411,37.8292331 137.948798,38.0775475 c 137.983078,38.4577806 137.927766,38.9415367 138.126037,39.2217479 c 138.240713,39.3838164 138.438596,39.433835 138.560544,39.6201964 c 138.747066,39.9052408 138.79591,40.4584527 138.951599,40.833933 c 139.061131,41.0980922 139.263918,41.2045958 139.429556,41.3702343 c 139.592672,41.5333495 139.626412,41.9345363 139.622169,42.4245876 c 139.61832,42.8690625 139.583226,43.386641 139.585076,43.8659947 c 139.586579,44.2554594 139.612471,44.6196912 139.699324,44.8989806 c 139.798518,45.2179557 140.058533,45.3368515 140.249382,45.5162438 c 140.468372,45.7220877 140.667698,45.9527405 140.883714,46.2064196 c 141.287987,46.6811791 141.363563,46.8042662 141.758733,46.8889452 c 142.153903,46.9736241 142.389571,46.8889452 142.5773,46.8889452 c 142.919018,46.8889452 143.062586,46.5447924 143.170055,46.4373229 c 143.5244,46.0829776 144.30261,46.0849327 144.456107,45.9576928 c 144.609603,45.8304528 144.538252,45.788898 144.680169,45.6469811 c 144.867391,45.4597592 145.3408,45.3492395 145.921047,45.2919329 c 146.416879,45.2429633 146.990723,45.23285 147.530673,45.2469365 c 148.234599,45.2653008 148.880917,45.3247952 149.221661,45.3929442 c 149.793527,45.5073172 149.592293,45.9406794 150.311344,46.0844896 c 150.590615,46.1403437 150.983508,46.1328596 151.385705,46.1339833 c 151.740461,46.1349745 152.102455,46.1426626 152.400101,46.2064196 c 153.035195,46.3424596 152.843259,46.6989649 153.303797,46.8553563 c 153.740109,47.0035208 154.126659,47.0607177 154.462767,47.0525307 c 154.848858,47.0431263 155.168389,46.9474451 155.420329,46.8042662 c 155.701153,46.6446723 155.856582,46.3018852 156.238444,46.2064196 c 156.63755,46.106643 157.062463,46.0773864 157.501193,46.0799003 c 157.962514,46.0825436 158.439111,46.1203133 158.917047,46.1481609 c 159.395158,46.1760187 159.874609,46.193947 160.341445,46.1568477 c 160.695624,46.1287011 161.042543,46.0688809 161.376107,45.9576928 c 161.556165,45.8976734 161.595501,45.6542159 161.855957,45.5162438 c 162.214357,45.3263872 162.76045,45.2045722 162.999127,45.1250133 c 163.301469,45.0242326 163.254211,44.8176843 163.542798,44.673391 c 163.811921,44.5388293 164.198751,44.5040294 164.537469,44.3628991 c 164.876187,44.2217687 165.245628,43.7958718 165.59596,43.4455407 c 165.694074,43.3474263 166.089427,42.951577 166.414031,42.5846337 c 166.738635,42.2176903 166.649398,41.6958381 166.82381,41.3702343 c 166.942587,41.1484928 167.231728,41.0717202 167.346164,40.833933 c 167.514797,40.4835286 167.470564,39.9584339 167.642542,39.6201964 c 167.752332,39.4042681 168.064756,39.3304646 168.164732,39.1262328 c 168.337576,38.7731432 168.382235,38.2192425 168.517395,37.9548375 c 168.663802,37.6684294 168.946353,37.5366951 169.011358,37.4066857 c 169.16199,37.1054214 169.196418,36.5836194 169.216001,36.0534301 c 169.23683,35.4895055 169.240866,34.9160924 169.350074,34.5884691 c 169.469774,34.2293693 169.81581,34.2906902 169.971058,33.7685027 c 170.126306,33.2463152 169.678613,32.7676354 169.40929,32.4983129 c 169.209305,32.2983273 169.341035,31.4437653 168.57385,30.8470659 c 168.319811,30.6494805 168.456555,30.0086667 168.164732,29.6053234 c 167.993437,29.3685691 167.724169,29.2856739 167.642542,29.0407946 c 167.455896,28.4808568 167.609275,28.4279622 167.346164,27.982304 c 167.124194,27.6063293 166.62278,27.1459759 166.287673,26.8554088 c 165.944635,26.5579643 165.711054,26.2484947 165.454331,26.0064531 c 165.144628,25.7144613 165.045061,25.484594 164.678766,25.2711034 c 164.34227,25.0749805 163.922398,25.0466617 163.542798,24.8773955 c 163.194478,24.7220773 163.083972,24.3552062 162.507536,24.2281877 c 161.9311,24.1011693 161.486505,24.2278936 161.023442,24.0729418 c 160.683926,23.9593313 160.613722,23.5600169 160.247214,23.4378478 c 159.788032,23.284787 159.142879,23.3079033 158.470122,23.3479115 c 157.944612,23.379163 157.40226,23.4207213 156.918545,23.3966683 c 156.540507,23.3778701 156.198284,23.318997 155.927908,23.1838091 c 155.706656,23.0731829 155.491061,22.6254308 154.912418,22.5487133 c 154.333776,22.4719957 153.815269,22.5408899 153.477743,22.3953734 c 153.162857,22.2596177 153.059763,21.8669612 152.90768,21.8289403 c 152.619885,21.7569915 152.344296,21.7089837 152.079711,21.6777671 c 151.723489,21.6357386 151.387213,21.6241467 151.067951,21.6255422 c 150.539712,21.6278512 150.05805,21.6657143 149.609675,21.6600961 c 149.164256,21.6545149 148.751688,21.606024 148.358942,21.4371413 c 147.975987,21.272469 147.898826,21.0950527 147.20865,20.4176183 c 146.996476,20.2093614 146.799175,20.0264466 146.620482,19.8623133 c 146.217909,19.4925413 145.909784,19.2180944 145.738826,18.9639568 c 145.491984,18.5970134 145.583581,18.1453902 145.368182,17.7925597 c 145.152783,17.4397292 145.165926,17.7984953 144.896094,16.9889995 c 144.810593,16.7324952 145.287599,16.6014434 145.368182,16.2988237 c 145.639372,15.2804067 145.643957,15.1850966 146.165982,14.6594246 c 146.360127,14.4639228 146.490629,13.4976765 146.72609,13.2622164 c 146.958739,13.0295674 147.222439,12.8230921 147.459269,12.5905543 c 147.691908,12.3621314 147.898621,12.1085601 148.024506,11.7803283 c 148.278545,11.1179515 147.982167,11.0791053 148.358942,10.9053089 c 148.735717,10.7315125 148.771309,11.160942 149.221661,11.3992716 c 149.504567,11.5489873 149.822314,11.5203687 150.197941,11.5545166 c 150.54641,11.5861956 151.129297,11.5778787 151.716397,11.5676281 c 152.606024,11.5520956 153.505322,11.5321234 153.61334,11.6401409 c 153.766508,11.7933085 153.971577,12.0864654 154.206096,12.2037249 c 154.698315,12.4498344 155.30673,12.4266608 155.920073,12.3855933 c 156.614755,12.3390797 157.315759,12.2696115 157.86142,12.5424421 c 157.984728,12.6040962 158.233373,12.9583493 158.341468,13.0123969 c 159.024284,13.3538048 159.278204,13.2615788 159.639685,13.4588758 c 159.942308,13.6240473 160.113559,13.8965364 160.486478,14.2710927 c 160.686801,14.4722946 160.886652,14.685981 161.1847,14.9840285 c 161.365783,15.1651113 161.695765,15.5370768 162.058957,15.9326287 c 162.26739,16.1596336 162.486761,16.3944067 162.695202,16.6053267 c 162.98153,16.8950608 163.247234,17.1397848 163.435633,17.2575338 c 164.000162,17.6103643 164.021282,17.325205 164.339387,17.4115026 c 164.810238,17.5392378 165.264047,17.4414623 165.657123,17.1984805 c 165.99567,16.9892063 166.289166,16.6722189 166.509697,16.2988237 c 166.701958,15.973295 166.74103,15.4880551 166.734104,14.9139927 c 166.728084,14.4149938 166.68731,13.8488822 166.682182,13.2622164 c 166.679036,12.9023848 166.695197,12.5006898 166.687695,12.1294325 c 166.677432,11.6215405 166.622882,11.170613 166.414031,10.9617617 c 166.128028,10.6757592 166.159076,10.7905709 165.962407,10.3972329 c 165.856212,10.1848444 165.924877,9.59242206 165.454331,9.16938297 c 164.983785,8.74634388 164.974481,8.36180797 164.452292,8.26613685 c 164.066307,8.19542002 163.685223,8.14909852 163.306397,8.12011209 c 162.791973,8.08075023 162.281711,8.07335435 161.768994,8.08024528 c 161.296026,8.08660196 160.820969,8.10511607 160.338626,8.12190977 c 159.851417,8.13887291 159.356775,8.15408076 158.849345,8.15323116 c 158.30471,8.15231927 157.5235,8.19792475 157.02874,8.07423488 c 156.478673,7.93671794 156.410757,7.53795352 155.927908,7.3770041 c 155.448624,7.21724255 154.659909,7.3770041 153.7078,7.27887313 c 153.241169,7.23077885 153.177424,6.91365194 152.90768,6.7284576 c 152.637688,6.54309314 152.23055,6.41794086 151.566924,6.47441966 c 151.12078,6.51238933 150.83207,6.57321218 150.371953,7.12362771 c 149.802771,7.36355241 149.279692,7.3710543 149.012431,7.50468461 c 148.879089,7.94916158 148.434746,8.09334816 147.953018,8.15070193 c 147.364023,8.22082663 146.719139,8.16114322 146.518475,8.36180796 c 146.394821,8.48546169 146.314956,8.65382226 146.165982,8.80279659 c 145.948678,9.02010087 145.440191,8.89079158 145.130308,9.04573276 c 144.88043,9.17067194 144.66026,9.44863746 144.456107,9.65279072 c 144.346843,9.76205412 143.798172,9.67084492 143.304335,10.0297316 c 142.810498,10.3886184 142.157865,11.0432512 142.08985,11.1355582 c 141.957628,11.3150011 141.643424,11.5108743 141.5439,11.7099224 c 141.102752,12.5922176 141.357822,12.4703809 141.236493,12.6980144 c 141.176405,12.8107505 141.061358,13.0475396 140.650022,13.4588759 c 140.642454,13.4664439 140.475326,13.8671445 140.479441,14.2710924 c 140.483232,14.6433395 140.496244,15.0055832 140.409424,15.1792235 c 140.115401,15.7672686 140.079292,15.5288997 139.848007,15.9935344 c 139.757941,16.1744701 139.58036,16.6889278 139.699324,16.9889992 c 139.885854,17.4594979 140.374476,17.9028886 140.409424,18.0776274 c 140.47668,18.413908 140.479334,18.662833 140.479441,19.1170398 c 140.479569,19.6665979 140.739019,19.9621383 141.05337,20.2266245 c 141.344686,20.4717288 141.683151,20.6901635 141.906073,21.0593792 c 142.369547,21.8270067 142.205884,22.3311986 142.369547,22.3953734 c 142.760603,22.5487133 143.328866,22.4966556 143.600649,22.768438 c 143.883821,23.05161 143.832385,23.9995403 144.093089,24.0729418 c 144.78673,24.2682369 144.91845,24.3326604 145.266426,24.4485295 c 145.614403,24.5643987 145.688025,24.9264879 146.165982,25.3320262 c 146.643938,25.7375646 146.888584,26.0064531 147.20865,26.3748412 c 147.473576,26.6797634 147.82896,26.7022626 148.178513,26.7926507 c 148.453283,26.863701 148.724449,26.9766989 148.945247,27.3017888 c 149.105281,27.5374137 149.59224,27.5897622 150.174708,27.5910411 c 150.625689,27.5920313 151.133927,27.5624062 151.592007,27.5635302 c 152.05168,27.5646582 152.460847,27.5967494 152.71097,27.7218108 c 153.119776,27.9262138 152.94786,27.9045245 153.190169,28.1468334 c 153.503682,28.460346 154.073659,28.4569839 154.654062,28.4357353 c 155.185456,28.416281 155.72559,28.3818333 156.085636,28.5618566 c 156.437966,28.7380211 156.369805,28.9284363 156.940166,29.2136168 c 157.074324,29.280696 157.497264,29.267697 157.908848,29.2925183 c 158.186808,29.3092812 158.459588,29.3432935 158.634741,29.4308696 c 159.069248,29.6481224 159.03255,29.8598897 159.195756,29.9142919 c 159.579797,30.0423054 160.386616,30.1669654 160.609554,30.3899032 c 160.99896,30.7793095 161.359188,31.2508628 161.647319,31.5165013 c 161.987582,31.8302031 162.21238,32.061759 162.332537,32.513503 c 162.403703,32.7810623 162.438162,33.1258629 162.438162,33.6107207 c 162.438162,34.0186619 162.452837,34.7374085 162.226661,34.9635845 c 162.099415,35.0908296 161.34304,35.1740848 161.184699,35.2532555 c 160.803153,35.4440287 160.436636,35.8990815 160.126547,36.2091704 c 159.841794,36.4939233 159.179177,37.2442251 158.849344,37.5271727 c 158.562323,37.7733947 157.94583,37.6887831 157.649859,37.8458129 c 157.288559,38.0375035 157.304334,38.2294232 157.028739,38.3672203 c 156.779615,38.4917827 156.349232,38.5325445 155.869915,38.5444681 c 155.482441,38.554107 155.062988,38.5449005 154.68146,38.5458838 c 154.13189,38.5473002 153.66101,38.5698593 153.477742,38.7003405 c 153.089172,38.9769904 153.019202,39.246292 152.448725,39.3350424 c 152.181539,39.3766091 151.804564,39.3785708 151.233648,39.3231334 c 150.418741,39.2440039 150.596106,38.8596605 150.197941,38.7003405 c 150.029005,38.632743 149.680384,38.5953741 149.321033,38.5608198 c 148.833433,38.5139334 148.326079,38.4722292 148.22107,38.3672203 c 147.962814,38.1089646 148.045898,38.0767299 147.786563,37.9548375 c 147.475527,37.8086441 146.977945,37.7125517 146.606729,37.4677074 c 146.218796,37.2118372 145.84714,36.8106482 145.465167,36.4029632 c 145.176871,36.0952613 144.882699,35.7838589 144.571216,35.5284429 c 144.330387,35.3309627 143.73117,35.1690929 143.304335,35.0360027 c 142.8775,34.9029125 143.0809,34.0577784 142.873936,33.8154271 c 142.600667,33.4954338 142.148359,33.3790233 141.446919,33.4138815 c 141.04774,33.4337188 140.567878,33.5025458 139.994415,33.6107207 c 139.62487,33.6804296 139.617843,34.0945714 138.951599,34.5884691 l 138.811517,35.1061011 z ";
// var path = point1 + point2 + point3 + point4 + point5 + point6 + point7 + point8 + point9;
var path = lettresANT + lettreS;

module.exports = path;
},{}],"/Users/Romain/Documents/Programmation/ants/AntColony/src/utilities.js":[function(require,module,exports){
'use strict';

var sqrt = Math.sqrt;
var pow = Math.pow;

function sign(x) {
	return x ? x < 0 ? -1 : 1 : 0;
}

function range(start, count) {
    return Array.apply(0, Array(count)).map(function (element, index) {
    	return index + start
    });
}

function distance(a, b){
	return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2));
}

function norm(v){
	return sqrt(pow(v.x, 2) + pow(v.y, 2));
}

module.exports = {
	sign: sign,
	range: range,
	distance: distance,
	norm: norm
}
},{}],"/Users/Romain/Documents/Programmation/ants/AntColony/src/vector.js":[function(require,module,exports){
'use strict'

function Vector(x, y) {
    this.x = x;                
    this.y = y;
}

Vector.prototype.norm = function(){
	return Math.sqrt(this.x * this.x + this.y * this.y);
}

Vector.prototype.normalize = function(){
	var norm = this.norm();
	this.x = this.x / norm;
	this.y = this.y / norm;
}



module.exports = Vector;
},{}],"/Users/Romain/Documents/Programmation/ants/AntColony/start.js":[function(require,module,exports){
'use strict';

var _antColony = require('./index.js');

var container = document.querySelector('.colony');

var options = {
	velocity: 0.001,
	nbAnts: 4000,
	weight: 10,
	repSize: 0.05,
	repSpeed: 0.002,
	nbStart: 300,
	nbRand: 300
	// obj par defaut
};

var antColony = _antColony(container, options);

window.addEventListener('click', function (){
	// options.velocity = 0.003;
	options.nbAnts = 20000;
	// options.weight = 10000000;
	// options.repSpeed = 0.01;
	// options.repSize = 0.1;

	// antColony.changeOptions(options);
	antColony.changeOptions(options);
});


},{"./index.js":"/Users/Romain/Documents/Programmation/ants/AntColony/index.js"}]},{},["/Users/Romain/Documents/Programmation/ants/AntColony/start.js"])
//# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJzb3VyY2VzIjpbIi91c3IvbG9jYWwvbGliL25vZGVfbW9kdWxlcy93YXRjaGlmeS9ub2RlX21vZHVsZXMvYnJvd3NlcmlmeS9ub2RlX21vZHVsZXMvYnJvd3Nlci1wYWNrL19wcmVsdWRlLmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9hbnRzL0FudENvbG9ueS9pbmRleC5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vYW50cy9BbnRDb2xvbnkvbm9kZV9tb2R1bGVzL2RlbGF1bmF5LXRyaWFuZ3VsYXRlL25vZGVfbW9kdWxlcy9pbmNyZW1lbnRhbC1jb252ZXgtaHVsbC9pY2guanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL2FudHMvQW50Q29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvaW5jcmVtZW50YWwtY29udmV4LWh1bGwvbm9kZV9tb2R1bGVzL3JvYnVzdC1vcmllbnRhdGlvbi9ub2RlX21vZHVsZXMvcm9idXN0LXNjYWxlL25vZGVfbW9kdWxlcy90d28tc3VtL3R3by1zdW0uanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL2FudHMvQW50Q29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvaW5jcmVtZW50YWwtY29udmV4LWh1bGwvbm9kZV9tb2R1bGVzL3JvYnVzdC1vcmllbnRhdGlvbi9ub2RlX21vZHVsZXMvcm9idXN0LXNjYWxlL3JvYnVzdC1zY2FsZS5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vYW50cy9BbnRDb2xvbnkvbm9kZV9tb2R1bGVzL2RlbGF1bmF5LXRyaWFuZ3VsYXRlL25vZGVfbW9kdWxlcy9pbmNyZW1lbnRhbC1jb252ZXgtaHVsbC9ub2RlX21vZHVsZXMvcm9idXN0LW9yaWVudGF0aW9uL25vZGVfbW9kdWxlcy9yb2J1c3Qtc3VidHJhY3Qvcm9idXN0LWRpZmYuanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL2FudHMvQW50Q29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvaW5jcmVtZW50YWwtY29udmV4LWh1bGwvbm9kZV9tb2R1bGVzL3JvYnVzdC1vcmllbnRhdGlvbi9ub2RlX21vZHVsZXMvcm9idXN0LXN1bS9yb2J1c3Qtc3VtLmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9hbnRzL0FudENvbG9ueS9ub2RlX21vZHVsZXMvZGVsYXVuYXktdHJpYW5ndWxhdGUvbm9kZV9tb2R1bGVzL2luY3JlbWVudGFsLWNvbnZleC1odWxsL25vZGVfbW9kdWxlcy9yb2J1c3Qtb3JpZW50YXRpb24vbm9kZV9tb2R1bGVzL3R3by1wcm9kdWN0L3R3by1wcm9kdWN0LmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9hbnRzL0FudENvbG9ueS9ub2RlX21vZHVsZXMvZGVsYXVuYXktdHJpYW5ndWxhdGUvbm9kZV9tb2R1bGVzL2luY3JlbWVudGFsLWNvbnZleC1odWxsL25vZGVfbW9kdWxlcy9yb2J1c3Qtb3JpZW50YXRpb24vb3JpZW50YXRpb24uanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL2FudHMvQW50Q29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvaW5jcmVtZW50YWwtY29udmV4LWh1bGwvbm9kZV9tb2R1bGVzL3NpbXBsaWNpYWwtY29tcGxleC9ub2RlX21vZHVsZXMvYml0LXR3aWRkbGUvdHdpZGRsZS5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vYW50cy9BbnRDb2xvbnkvbm9kZV9tb2R1bGVzL2RlbGF1bmF5LXRyaWFuZ3VsYXRlL25vZGVfbW9kdWxlcy9pbmNyZW1lbnRhbC1jb252ZXgtaHVsbC9ub2RlX21vZHVsZXMvc2ltcGxpY2lhbC1jb21wbGV4L25vZGVfbW9kdWxlcy91bmlvbi1maW5kL2luZGV4LmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9hbnRzL0FudENvbG9ueS9ub2RlX21vZHVsZXMvZGVsYXVuYXktdHJpYW5ndWxhdGUvbm9kZV9tb2R1bGVzL2luY3JlbWVudGFsLWNvbnZleC1odWxsL25vZGVfbW9kdWxlcy9zaW1wbGljaWFsLWNvbXBsZXgvdG9wb2xvZ3kuanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL2FudHMvQW50Q29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvdW5pcS91bmlxLmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9hbnRzL0FudENvbG9ueS9ub2RlX21vZHVsZXMvZGVsYXVuYXktdHJpYW5ndWxhdGUvdHJpYW5ndWxhdGUuanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL2FudHMvQW50Q29sb255L25vZGVfbW9kdWxlcy9wYXJzZS1zdmctcGF0aC9pbmRleC5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vYW50cy9BbnRDb2xvbnkvc3JjL2FudC5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vYW50cy9BbnRDb2xvbnkvc3JjL2FudHNHcm91cC5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vYW50cy9BbnRDb2xvbnkvc3JjL2NyZWF0ZUVkZ2VzLmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9hbnRzL0FudENvbG9ueS9zcmMvZWRnZS5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vYW50cy9BbnRDb2xvbnkvc3JjL2luaXRpYWxpemVQb2ludHMuanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL2FudHMvQW50Q29sb255L3NyYy9tb3VzZS5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vYW50cy9BbnRDb2xvbnkvc3JjL3BvaW50LmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9hbnRzL0FudENvbG9ueS9zcmMvcmVuZGVyaW5nLmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9hbnRzL0FudENvbG9ueS9zcmMvc3ZnUGF0aC5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vYW50cy9BbnRDb2xvbnkvc3JjL3V0aWxpdGllcy5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vYW50cy9BbnRDb2xvbnkvc3JjL3ZlY3Rvci5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vYW50cy9BbnRDb2xvbnkvc3RhcnQuanMiXSwibmFtZXMiOltdLCJtYXBwaW5ncyI6IkFBQUE7QUNBQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUN2Q0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUM3YkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNoQkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNqREE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQzNKQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDM0pBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNoQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDN0xBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQzVNQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUN6REE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDdFZBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ3pEQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDOUpBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDdkRBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDbk5BO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNuQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNyRkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDOUVBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ3hLQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ25CQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNUQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDNU9BO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDaEJBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDNUJBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDbkJBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBIiwiZmlsZSI6ImdlbmVyYXRlZC5qcyIsInNvdXJjZVJvb3QiOiIiLCJzb3VyY2VzQ29udGVudCI6WyIoZnVuY3Rpb24gZSh0LG4scil7ZnVuY3Rpb24gcyhvLHUpe2lmKCFuW29dKXtpZighdFtvXSl7dmFyIGE9dHlwZW9mIHJlcXVpcmU9PVwiZnVuY3Rpb25cIiYmcmVxdWlyZTtpZighdSYmYSlyZXR1cm4gYShvLCEwKTtpZihpKXJldHVybiBpKG8sITApO3ZhciBmPW5ldyBFcnJvcihcIkNhbm5vdCBmaW5kIG1vZHVsZSAnXCIrbytcIidcIik7dGhyb3cgZi5jb2RlPVwiTU9EVUxFX05PVF9GT1VORFwiLGZ9dmFyIGw9bltvXT17ZXhwb3J0czp7fX07dFtvXVswXS5jYWxsKGwuZXhwb3J0cyxmdW5jdGlvbihlKXt2YXIgbj10W29dWzFdW2VdO3JldHVybiBzKG4/bjplKX0sbCxsLmV4cG9ydHMsZSx0LG4scil9cmV0dXJuIG5bb10uZXhwb3J0c312YXIgaT10eXBlb2YgcmVxdWlyZT09XCJmdW5jdGlvblwiJiZyZXF1aXJlO2Zvcih2YXIgbz0wO288ci5sZW5ndGg7bysrKXMocltvXSk7cmV0dXJuIHN9KSIsIid1c2Ugc3RyaWN0JztcblxudmFyIGluaXRSZW5kZXJpbmcgPSByZXF1aXJlKCcuL3NyYy9yZW5kZXJpbmcuanMnKTtcbnZhciBpbml0aWFsaXplUG9pbnRzID0gcmVxdWlyZSgnLi9zcmMvaW5pdGlhbGl6ZVBvaW50cy5qcycpO1xudmFyIGNyZWF0ZUVkZ2VzID0gcmVxdWlyZSgnLi9zcmMvY3JlYXRlRWRnZXMuanMnKTtcbi8vIHZhciBpbml0QW50cyA9IHJlcXVpcmUoJy4vc3JjL2luaXRpYWxpemVBbnRzJyk7XG5cbm1vZHVsZS5leHBvcnRzID0gZnVuY3Rpb24gaW5pdChjb250YWluZXJFbGVtZW50LCBvcHRpb25zKXtcblxuXHR2YXIgcmVuZGVyLCBwb2ludHNJbmZvcywgZWRnZXMsIHBvcHVsYXRpb24sIHBvaW50c01hcDtcblxuXG5cdGZ1bmN0aW9uIF9pbml0KGNvbnRhaW5lckVsZW1lbnQsIG9wdGlvbnMpe1xuXHRcdHBvaW50c0luZm9zID0gaW5pdGlhbGl6ZVBvaW50cyhvcHRpb25zLm5iU3RhcnQsIG9wdGlvbnMubmJSYW5kKTtcblx0XHRlZGdlcyA9IGNyZWF0ZUVkZ2VzKHBvaW50c0luZm9zLnBvaW50cyk7XG5cdFx0Ly8gcG9wdWxhdGlvbiA9IG9wdGlvbnMubmJBbnRzO1xuXHRcdC8vIHBvcHVsYXRpb24gPSBpbml0QW50cyhjb250YWluZXJFbGVtZW50LCBwb2ludHNJbmZvcywgb3B0aW9ucyk7XG5cdFx0cG9pbnRzTWFwID0ge1xuXHRcdFx0cG9pbnRzSW5mb3M6IHBvaW50c0luZm9zLFxuXHRcdFx0ZWRnZXM6IGVkZ2VzXG5cdFx0XHQvLyBwb3B1bGF0aW9uOiBwb3B1bGF0aW9uXG5cdFx0fTtcblx0XHRyZW5kZXIgPSBpbml0UmVuZGVyaW5nKGNvbnRhaW5lckVsZW1lbnQsIHBvaW50c01hcCwgb3B0aW9ucyk7XG5cdH1cblxuXHRfaW5pdChjb250YWluZXJFbGVtZW50LCBvcHRpb25zKTtcblxuXHRyZXR1cm4ge1xuXHRcdHRvZ2dsZVBsYXlQYXVzZTogZnVuY3Rpb24oKXsgcmVuZGVyLnRvZ2dsZVBsYXlQYXVzZSgpIH0sXG5cdFx0Y2hhbmdlT3B0aW9uczogZnVuY3Rpb24ob3B0cyl7XG5cdFx0XHRyZW5kZXIubW9kaWZ5QW50cyhvcHRzKTtcblx0XHR9LFxuXHRcdHJlc2V0OiBmdW5jdGlvbihvcHRzKXtcblx0XHRcdHJlbmRlci5yZXNldCgpO1xuXG5cdFx0XHRcdC8vIHJlc2V0IGVsZW1lbnRzXG5cdFx0XHRfaW5pdChjb250YWluZXJFbGVtZW50LCBvcHRzKTtcblx0XHR9XG5cdH07XG59OyIsIlwidXNlIHN0cmljdFwiXG5cbi8vSGlnaCBsZXZlbCBpZGVhOlxuLy8gMS4gVXNlIENsYXJrc29uJ3MgaW5jcmVtZW50YWwgY29uc3RydWN0aW9uIHRvIGZpbmQgY29udmV4IGh1bGxcbi8vIDIuIFBvaW50IGxvY2F0aW9uIGluIHRyaWFuZ3VsYXRpb24gYnkganVtcCBhbmQgd2Fsa1xuXG5tb2R1bGUuZXhwb3J0cyA9IGluY3JlbWVudGFsQ29udmV4SHVsbFxuXG52YXIgb3JpZW50ID0gcmVxdWlyZShcInJvYnVzdC1vcmllbnRhdGlvblwiKVxudmFyIGNvbXBhcmVDZWxsID0gcmVxdWlyZShcInNpbXBsaWNpYWwtY29tcGxleFwiKS5jb21wYXJlQ2VsbHNcblxuZnVuY3Rpb24gY29tcGFyZUludChhLCBiKSB7XG4gIHJldHVybiBhIC0gYlxufVxuXG5mdW5jdGlvbiBTaW1wbGV4KHZlcnRpY2VzLCBhZGphY2VudCwgYm91bmRhcnkpIHtcbiAgdGhpcy52ZXJ0aWNlcyA9IHZlcnRpY2VzXG4gIHRoaXMuYWRqYWNlbnQgPSBhZGphY2VudFxuICB0aGlzLmJvdW5kYXJ5ID0gYm91bmRhcnlcbiAgdGhpcy5sYXN0VmlzaXRlZCA9IC0xXG59XG5cblNpbXBsZXgucHJvdG90eXBlLmZsaXAgPSBmdW5jdGlvbigpIHtcbiAgdmFyIHQgPSB0aGlzLnZlcnRpY2VzWzBdXG4gIHRoaXMudmVydGljZXNbMF0gPSB0aGlzLnZlcnRpY2VzWzFdXG4gIHRoaXMudmVydGljZXNbMV0gPSB0XG4gIHZhciB1ID0gdGhpcy5hZGphY2VudFswXVxuICB0aGlzLmFkamFjZW50WzBdID0gdGhpcy5hZGphY2VudFsxXVxuICB0aGlzLmFkamFjZW50WzFdID0gdVxufVxuXG5mdW5jdGlvbiBHbHVlRmFjZXQodmVydGljZXMsIGNlbGwsIGluZGV4KSB7XG4gIHRoaXMudmVydGljZXMgPSB2ZXJ0aWNlc1xuICB0aGlzLmNlbGwgPSBjZWxsXG4gIHRoaXMuaW5kZXggPSBpbmRleFxufVxuXG5mdW5jdGlvbiBjb21wYXJlR2x1ZShhLCBiKSB7XG4gIHJldHVybiBjb21wYXJlQ2VsbChhLnZlcnRpY2VzLCBiLnZlcnRpY2VzKVxufVxuXG5mdW5jdGlvbiBiYWtlT3JpZW50KGQpIHtcbiAgdmFyIGNvZGUgPSBbXCJmdW5jdGlvbiBvcmllbnQoKXt2YXIgdHVwbGU9dGhpcy50dXBsZTtyZXR1cm4gdGVzdChcIl1cbiAgZm9yKHZhciBpPTA7IGk8PWQ7ICsraSkge1xuICAgIGlmKGkgPiAwKSB7XG4gICAgICBjb2RlLnB1c2goXCIsXCIpXG4gICAgfVxuICAgIGNvZGUucHVzaChcInR1cGxlW1wiLCBpLCBcIl1cIilcbiAgfVxuICBjb2RlLnB1c2goXCIpfXJldHVybiBvcmllbnRcIilcbiAgdmFyIHByb2MgPSBuZXcgRnVuY3Rpb24oXCJ0ZXN0XCIsIGNvZGUuam9pbihcIlwiKSlcbiAgdmFyIHRlc3QgPSBvcmllbnRbZCsxXVxuICBpZighdGVzdCkge1xuICAgIHRlc3QgPSBvcmllbnRcbiAgfVxuICByZXR1cm4gcHJvYyh0ZXN0KVxufVxuXG52YXIgQkFLRUQgPSBbXVxuXG5mdW5jdGlvbiBUcmlhbmd1bGF0aW9uKGRpbWVuc2lvbiwgdmVydGljZXMsIHNpbXBsaWNlcykge1xuICB0aGlzLmRpbWVuc2lvbiA9IGRpbWVuc2lvblxuICB0aGlzLnZlcnRpY2VzID0gdmVydGljZXNcbiAgdGhpcy5zaW1wbGljZXMgPSBzaW1wbGljZXNcbiAgdGhpcy5pbnRlcmlvciA9IHNpbXBsaWNlcy5maWx0ZXIoZnVuY3Rpb24oYykge1xuICAgIHJldHVybiAhYy5ib3VuZGFyeVxuICB9KVxuXG4gIHRoaXMudHVwbGUgPSBuZXcgQXJyYXkoZGltZW5zaW9uKzEpXG4gIGZvcih2YXIgaT0wOyBpPD1kaW1lbnNpb247ICsraSkge1xuICAgIHRoaXMudHVwbGVbaV0gPSB0aGlzLnZlcnRpY2VzW2ldXG4gIH1cblxuICB2YXIgbyA9IEJBS0VEW2RpbWVuc2lvbl1cbiAgaWYoIW8pIHtcbiAgICBvID0gQkFLRURbZGltZW5zaW9uXSA9IGJha2VPcmllbnQoZGltZW5zaW9uKVxuICB9XG4gIHRoaXMub3JpZW50ID0gb1xufVxuXG52YXIgcHJvdG8gPSBUcmlhbmd1bGF0aW9uLnByb3RvdHlwZVxuXG4vL0RlZ2VuZXJhdGUgc2l0dWF0aW9uIHdoZXJlIHdlIGFyZSBvbiBib3VuZGFyeSwgYnV0IGNvcGxhbmFyIHRvIGZhY2VcbnByb3RvLmhhbmRsZUJvdW5kYXJ5RGVnZW5lcmFjeSA9IGZ1bmN0aW9uKGNlbGwsIHBvaW50KSB7XG4gIHZhciBkID0gdGhpcy5kaW1lbnNpb25cbiAgdmFyIG4gPSB0aGlzLnZlcnRpY2VzLmxlbmd0aCAtIDFcbiAgdmFyIHR1cGxlID0gdGhpcy50dXBsZVxuICB2YXIgdmVydHMgPSB0aGlzLnZlcnRpY2VzXG5cbiAgLy9EdW1iIHNvbHV0aW9uOiBKdXN0IGRvIGRmcyBmcm9tIGJvdW5kYXJ5IGNlbGwgdW50aWwgd2UgZmluZCBhbnkgcGVhaywgb3IgdGVybWluYXRlXG4gIHZhciB0b1Zpc2l0ID0gWyBjZWxsIF1cbiAgY2VsbC5sYXN0VmlzaXRlZCA9IC1uXG4gIHdoaWxlKHRvVmlzaXQubGVuZ3RoID4gMCkge1xuICAgIGNlbGwgPSB0b1Zpc2l0LnBvcCgpXG4gICAgdmFyIGNlbGxWZXJ0cyA9IGNlbGwudmVydGljZXNcbiAgICB2YXIgY2VsbEFkaiA9IGNlbGwuYWRqYWNlbnRcbiAgICBmb3IodmFyIGk9MDsgaTw9ZDsgKytpKSB7XG4gICAgICB2YXIgbmVpZ2hib3IgPSBjZWxsQWRqW2ldXG4gICAgICBpZighbmVpZ2hib3IuYm91bmRhcnkgfHwgbmVpZ2hib3IubGFzdFZpc2l0ZWQgPD0gLW4pIHtcbiAgICAgICAgY29udGludWVcbiAgICAgIH1cbiAgICAgIHZhciBudiA9IG5laWdoYm9yLnZlcnRpY2VzXG4gICAgICBmb3IodmFyIGo9MDsgajw9ZDsgKytqKSB7XG4gICAgICAgIHZhciB2diA9IG52W2pdXG4gICAgICAgIGlmKHZ2IDwgMCkge1xuICAgICAgICAgIHR1cGxlW2pdID0gcG9pbnRcbiAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICB0dXBsZVtqXSA9IHZlcnRzW3Z2XVxuICAgICAgICB9XG4gICAgICB9XG4gICAgICB2YXIgbyA9IHRoaXMub3JpZW50KClcbiAgICAgIGlmKG8gPiAwKSB7XG4gICAgICAgIHJldHVybiBuZWlnaGJvclxuICAgICAgfVxuICAgICAgbmVpZ2hib3IubGFzdFZpc2l0ZWQgPSAtblxuICAgICAgaWYobyA9PT0gMCkge1xuICAgICAgICB0b1Zpc2l0LnB1c2gobmVpZ2hib3IpXG4gICAgICB9XG4gICAgfVxuICB9XG4gIHJldHVybiBudWxsXG59XG5cbnByb3RvLndhbGsgPSBmdW5jdGlvbihwb2ludCwgcmFuZG9tKSB7XG4gIC8vQWxpYXMgbG9jYWwgcHJvcGVydGllc1xuICB2YXIgbiA9IHRoaXMudmVydGljZXMubGVuZ3RoIC0gMVxuICB2YXIgZCA9IHRoaXMuZGltZW5zaW9uXG4gIHZhciB2ZXJ0cyA9IHRoaXMudmVydGljZXNcbiAgdmFyIHR1cGxlID0gdGhpcy50dXBsZVxuXG4gIC8vQ29tcHV0ZSBpbml0aWFsIGp1bXAgY2VsbFxuICB2YXIgaW5pdEluZGV4ID0gcmFuZG9tID8gKHRoaXMuaW50ZXJpb3IubGVuZ3RoICogTWF0aC5yYW5kb20oKSl8MCA6ICh0aGlzLmludGVyaW9yLmxlbmd0aC0xKVxuICB2YXIgY2VsbCA9IHRoaXMuaW50ZXJpb3JbIGluaXRJbmRleCBdXG5cbiAgLy9TdGFydCB3YWxraW5nXG5vdXRlckxvb3A6XG4gIHdoaWxlKCFjZWxsLmJvdW5kYXJ5KSB7XG4gICAgdmFyIGNlbGxWZXJ0cyA9IGNlbGwudmVydGljZXNcbiAgICB2YXIgY2VsbEFkaiA9IGNlbGwuYWRqYWNlbnRcblxuICAgIGZvcih2YXIgaT0wOyBpPD1kOyArK2kpIHtcbiAgICAgIHR1cGxlW2ldID0gdmVydHNbY2VsbFZlcnRzW2ldXVxuICAgIH1cbiAgICBjZWxsLmxhc3RWaXNpdGVkID0gblxuXG4gICAgLy9GaW5kIGZhcnRoZXN0IGFkamFjZW50IGNlbGxcbiAgICBmb3IodmFyIGk9MDsgaTw9ZDsgKytpKSB7XG4gICAgICB2YXIgbmVpZ2hib3IgPSBjZWxsQWRqW2ldXG4gICAgICBpZihuZWlnaGJvci5sYXN0VmlzaXRlZCA+PSBuKSB7XG4gICAgICAgIGNvbnRpbnVlXG4gICAgICB9XG4gICAgICB2YXIgcHJldiA9IHR1cGxlW2ldXG4gICAgICB0dXBsZVtpXSA9IHBvaW50XG4gICAgICB2YXIgbyA9IHRoaXMub3JpZW50KClcbiAgICAgIHR1cGxlW2ldID0gcHJldlxuICAgICAgaWYobyA8IDApIHtcbiAgICAgICAgY2VsbCA9IG5laWdoYm9yXG4gICAgICAgIGNvbnRpbnVlIG91dGVyTG9vcFxuICAgICAgfSBlbHNlIHtcbiAgICAgICAgaWYoIW5laWdoYm9yLmJvdW5kYXJ5KSB7XG4gICAgICAgICAgbmVpZ2hib3IubGFzdFZpc2l0ZWQgPSBuXG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgbmVpZ2hib3IubGFzdFZpc2l0ZWQgPSAtblxuICAgICAgICB9XG4gICAgICB9XG4gICAgfVxuICAgIHJldHVyblxuICB9XG5cbiAgcmV0dXJuIGNlbGxcbn1cblxucHJvdG8uYWRkUGVha3MgPSBmdW5jdGlvbihwb2ludCwgY2VsbCkge1xuICB2YXIgbiA9IHRoaXMudmVydGljZXMubGVuZ3RoIC0gMVxuICB2YXIgZCA9IHRoaXMuZGltZW5zaW9uXG4gIHZhciB2ZXJ0cyA9IHRoaXMudmVydGljZXNcbiAgdmFyIHR1cGxlID0gdGhpcy50dXBsZVxuICB2YXIgaW50ZXJpb3IgPSB0aGlzLmludGVyaW9yXG4gIHZhciBzaW1wbGljZXMgPSB0aGlzLnNpbXBsaWNlc1xuXG4gIC8vV2Fsa2luZyBmaW5pc2hlZCBhdCBib3VuZGFyeSwgdGltZSB0byBhZGQgcGVha3NcbiAgdmFyIHRvdmlzaXQgPSBbIGNlbGwgXVxuXG4gIC8vU3RyZXRjaCBpbml0aWFsIGJvdW5kYXJ5IGNlbGwgaW50byBhIHBlYWtcbiAgY2VsbC5sYXN0VmlzaXRlZCA9IG5cbiAgY2VsbC52ZXJ0aWNlc1tjZWxsLnZlcnRpY2VzLmluZGV4T2YoLTEpXSA9IG5cbiAgY2VsbC5ib3VuZGFyeSA9IGZhbHNlXG4gIGludGVyaW9yLnB1c2goY2VsbClcblxuICAvL1JlY29yZCBhIGxpc3Qgb2YgYWxsIG5ldyBib3VuZGFyaWVzIGNyZWF0ZWQgYnkgYWRkZWQgcGVha3Mgc28gd2UgY2FuIGdsdWUgdGhlbSB0b2dldGhlciB3aGVuIHdlIGFyZSBhbGwgZG9uZVxuICB2YXIgZ2x1ZUZhY2V0cyA9IFtdXG5cbiAgLy9EbyBhIHRyYXZlcnNhbCBvZiB0aGUgYm91bmRhcnkgd2Fsa2luZyBvdXR3YXJkIGZyb20gc3RhcnRpbmcgcGVha1xuICB3aGlsZSh0b3Zpc2l0Lmxlbmd0aCA+IDApIHtcbiAgICAvL1BvcCBvZmYgcGVhayBhbmQgd2FsayBvdmVyIGFkamFjZW50IGNlbGxzXG4gICAgdmFyIGNlbGwgPSB0b3Zpc2l0LnBvcCgpXG4gICAgdmFyIGNlbGxWZXJ0cyA9IGNlbGwudmVydGljZXNcbiAgICB2YXIgY2VsbEFkaiA9IGNlbGwuYWRqYWNlbnRcbiAgICB2YXIgaW5kZXhPZk4gPSBjZWxsVmVydHMuaW5kZXhPZihuKVxuICAgIGlmKGluZGV4T2ZOIDwgMCkge1xuICAgICAgY29udGludWVcbiAgICB9XG5cbiAgICBmb3IodmFyIGk9MDsgaTw9ZDsgKytpKSB7XG4gICAgICBpZihpID09PSBpbmRleE9mTikge1xuICAgICAgICBjb250aW51ZVxuICAgICAgfVxuXG4gICAgICAvL0ZvciBlYWNoIGJvdW5kYXJ5IG5laWdoYm9yIG9mIHRoZSBjZWxsXG4gICAgICB2YXIgbmVpZ2hib3IgPSBjZWxsQWRqW2ldXG4gICAgICBpZighbmVpZ2hib3IuYm91bmRhcnkgfHwgbmVpZ2hib3IubGFzdFZpc2l0ZWQgPj0gbikge1xuICAgICAgICBjb250aW51ZVxuICAgICAgfVxuXG4gICAgICB2YXIgbnYgPSBuZWlnaGJvci52ZXJ0aWNlc1xuXG4gICAgICAvL1Rlc3QgaWYgbmVpZ2hib3IgaXMgYSBwZWFrXG4gICAgICBpZihuZWlnaGJvci5sYXN0VmlzaXRlZCAhPT0gLW4pIHsgICAgICBcbiAgICAgICAgLy9Db21wdXRlIG9yaWVudGF0aW9uIG9mIHAgcmVsYXRpdmUgdG8gZWFjaCBib3VuZGFyeSBwZWFrXG4gICAgICAgIHZhciBpbmRleE9mTmVnMSA9IDBcbiAgICAgICAgZm9yKHZhciBqPTA7IGo8PWQ7ICsraikge1xuICAgICAgICAgIGlmKG52W2pdIDwgMCkge1xuICAgICAgICAgICAgaW5kZXhPZk5lZzEgPSBqXG4gICAgICAgICAgICB0dXBsZVtqXSA9IHBvaW50XG4gICAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgIHR1cGxlW2pdID0gdmVydHNbbnZbal1dXG4gICAgICAgICAgfVxuICAgICAgICB9XG4gICAgICAgIHZhciBvID0gdGhpcy5vcmllbnQoKVxuXG4gICAgICAgIC8vVGVzdCBpZiBuZWlnaGJvciBjZWxsIGlzIGFsc28gYSBwZWFrXG4gICAgICAgIGlmKG8gPiAwKSB7XG4gICAgICAgICAgbnZbaW5kZXhPZk5lZzFdID0gblxuICAgICAgICAgIG5laWdoYm9yLmJvdW5kYXJ5ID0gZmFsc2VcbiAgICAgICAgICBpbnRlcmlvci5wdXNoKG5laWdoYm9yKVxuICAgICAgICAgIHRvdmlzaXQucHVzaChuZWlnaGJvcilcbiAgICAgICAgICBuZWlnaGJvci5sYXN0VmlzaXRlZCA9IG5cbiAgICAgICAgICBjb250aW51ZVxuICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgIG5laWdoYm9yLmxhc3RWaXNpdGVkID0gLW5cbiAgICAgICAgfVxuICAgICAgfVxuXG4gICAgICB2YXIgbmEgPSBuZWlnaGJvci5hZGphY2VudFxuXG4gICAgICAvL090aGVyd2lzZSwgcmVwbGFjZSBuZWlnaGJvciB3aXRoIG5ldyBmYWNlXG4gICAgICB2YXIgdnZlcnRzID0gY2VsbFZlcnRzLnNsaWNlKClcbiAgICAgIHZhciB2YWRqID0gY2VsbEFkai5zbGljZSgpXG4gICAgICB2YXIgbmNlbGwgPSBuZXcgU2ltcGxleCh2dmVydHMsIHZhZGosIHRydWUpXG4gICAgICBzaW1wbGljZXMucHVzaChuY2VsbClcblxuICAgICAgLy9Db25uZWN0IHRvIG5laWdoYm9yXG4gICAgICB2YXIgb3Bwb3NpdGUgPSBuYS5pbmRleE9mKGNlbGwpXG4gICAgICBpZihvcHBvc2l0ZSA8IDApIHtcbiAgICAgICAgY29udGludWVcbiAgICAgIH1cbiAgICAgIG5hW29wcG9zaXRlXSA9IG5jZWxsXG4gICAgICB2YWRqW2luZGV4T2ZOXSA9IG5laWdoYm9yXG5cbiAgICAgIC8vQ29ubmVjdCB0byBjZWxsXG4gICAgICB2dmVydHNbaV0gPSAtMVxuICAgICAgdmFkaltpXSA9IGNlbGxcbiAgICAgIGNlbGxBZGpbaV0gPSBuY2VsbFxuXG4gICAgICAvL0ZsaXAgZmFjZXRcbiAgICAgIG5jZWxsLmZsaXAoKVxuXG4gICAgICAvL0FkZCB0byBnbHVlIGxpc3RcbiAgICAgIGZvcih2YXIgaj0wOyBqPD1kOyArK2opIHtcbiAgICAgICAgdmFyIHV1ID0gdnZlcnRzW2pdXG4gICAgICAgIGlmKHV1IDwgMCB8fCB1dSA9PT0gbikge1xuICAgICAgICAgIGNvbnRpbnVlXG4gICAgICAgIH1cbiAgICAgICAgdmFyIG5mYWNlID0gbmV3IEFycmF5KGQtMSlcbiAgICAgICAgdmFyIG5wdHIgPSAwXG4gICAgICAgIGZvcih2YXIgaz0wOyBrPD1kOyArK2spIHtcbiAgICAgICAgICB2YXIgdnYgPSB2dmVydHNba11cbiAgICAgICAgICBpZih2diA8IDAgfHwgayA9PT0gaikge1xuICAgICAgICAgICAgY29udGludWVcbiAgICAgICAgICB9XG4gICAgICAgICAgbmZhY2VbbnB0cisrXSA9IHZ2XG4gICAgICAgIH1cbiAgICAgICAgZ2x1ZUZhY2V0cy5wdXNoKG5ldyBHbHVlRmFjZXQobmZhY2UsIG5jZWxsLCBqKSlcbiAgICAgIH1cbiAgICB9XG4gIH1cblxuICAvL0dsdWUgYm91bmRhcnkgZmFjZXRzIHRvZ2V0aGVyXG4gIGdsdWVGYWNldHMuc29ydChjb21wYXJlR2x1ZSlcblxuICBmb3IodmFyIGk9MDsgaSsxPGdsdWVGYWNldHMubGVuZ3RoOyBpKz0yKSB7XG4gICAgdmFyIGEgPSBnbHVlRmFjZXRzW2ldXG4gICAgdmFyIGIgPSBnbHVlRmFjZXRzW2krMV1cbiAgICB2YXIgYWkgPSBhLmluZGV4XG4gICAgdmFyIGJpID0gYi5pbmRleFxuICAgIGlmKGFpIDwgMCB8fCBiaSA8IDApIHtcbiAgICAgIGNvbnRpbnVlXG4gICAgfVxuICAgIGEuY2VsbC5hZGphY2VudFthLmluZGV4XSA9IGIuY2VsbFxuICAgIGIuY2VsbC5hZGphY2VudFtiLmluZGV4XSA9IGEuY2VsbFxuICB9XG59XG5cbnByb3RvLmluc2VydCA9IGZ1bmN0aW9uKHBvaW50LCByYW5kb20pIHtcbiAgLy9BZGQgcG9pbnRcbiAgdmFyIHZlcnRzID0gdGhpcy52ZXJ0aWNlc1xuICB2ZXJ0cy5wdXNoKHBvaW50KVxuXG4gIHZhciBjZWxsID0gdGhpcy53YWxrKHBvaW50LCByYW5kb20pXG4gIGlmKCFjZWxsKSB7XG4gICAgcmV0dXJuXG4gIH1cblxuICAvL0FsaWFzIGxvY2FsIHByb3BlcnRpZXNcbiAgdmFyIGQgPSB0aGlzLmRpbWVuc2lvblxuICB2YXIgdHVwbGUgPSB0aGlzLnR1cGxlXG5cbiAgLy9EZWdlbmVyYXRlIGNhc2U6IElmIHBvaW50IGlzIGNvcGxhbmFyIHRvIGNlbGwsIHRoZW4gd2FsayB1bnRpbCB3ZSBmaW5kIGEgbm9uLWRlZ2VuZXJhdGUgYm91bmRhcnlcbiAgZm9yKHZhciBpPTA7IGk8PWQ7ICsraSkge1xuICAgIHZhciB2diA9IGNlbGwudmVydGljZXNbaV1cbiAgICBpZih2diA8IDApIHtcbiAgICAgIHR1cGxlW2ldID0gcG9pbnRcbiAgICB9IGVsc2Uge1xuICAgICAgdHVwbGVbaV0gPSB2ZXJ0c1t2dl1cbiAgICB9XG4gIH1cbiAgdmFyIG8gPSB0aGlzLm9yaWVudCh0dXBsZSlcbiAgaWYobyA8IDApIHtcbiAgICByZXR1cm5cbiAgfSBlbHNlIGlmKG8gPT09IDApIHtcbiAgICBjZWxsID0gdGhpcy5oYW5kbGVCb3VuZGFyeURlZ2VuZXJhY3koY2VsbCwgcG9pbnQpXG4gICAgaWYoIWNlbGwpIHtcbiAgICAgIHJldHVyblxuICAgIH1cbiAgfVxuXG4gIC8vQWRkIHBlYWtzXG4gIHRoaXMuYWRkUGVha3MocG9pbnQsIGNlbGwpXG59XG5cbi8vRXh0cmFjdCBhbGwgYm91bmRhcnkgY2VsbHNcbnByb3RvLmJvdW5kYXJ5ID0gZnVuY3Rpb24oKSB7XG4gIHZhciBkID0gdGhpcy5kaW1lbnNpb25cbiAgdmFyIGJvdW5kYXJ5ID0gW11cbiAgdmFyIGNlbGxzID0gdGhpcy5zaW1wbGljZXNcbiAgdmFyIG5jID0gY2VsbHMubGVuZ3RoXG4gIGZvcih2YXIgaT0wOyBpPG5jOyArK2kpIHtcbiAgICB2YXIgYyA9IGNlbGxzW2ldXG4gICAgaWYoYy5ib3VuZGFyeSkge1xuICAgICAgdmFyIGJjZWxsID0gbmV3IEFycmF5KGQpXG4gICAgICB2YXIgY3YgPSBjLnZlcnRpY2VzXG4gICAgICB2YXIgcHRyID0gMFxuICAgICAgdmFyIHBhcml0eSA9IDBcbiAgICAgIGZvcih2YXIgaj0wOyBqPD1kOyArK2opIHtcbiAgICAgICAgaWYoY3Zbal0gPj0gMCkge1xuICAgICAgICAgIGJjZWxsW3B0cisrXSA9IGN2W2pdXG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgcGFyaXR5ID0gaiYxXG4gICAgICAgIH1cbiAgICAgIH1cbiAgICAgIGlmKHBhcml0eSA9PT0gKGQmMSkpIHtcbiAgICAgICAgdmFyIHQgPSBiY2VsbFswXVxuICAgICAgICBiY2VsbFswXSA9IGJjZWxsWzFdXG4gICAgICAgIGJjZWxsWzFdID0gdFxuICAgICAgfVxuICAgICAgYm91bmRhcnkucHVzaChiY2VsbClcbiAgICB9XG4gIH1cbiAgcmV0dXJuIGJvdW5kYXJ5XG59XG5cbmZ1bmN0aW9uIGluY3JlbWVudGFsQ29udmV4SHVsbChwb2ludHMsIHJhbmRvbVNlYXJjaCkge1xuICB2YXIgbiA9IHBvaW50cy5sZW5ndGhcbiAgaWYobiA9PT0gMCkge1xuICAgIHRocm93IG5ldyBFcnJvcihcIk11c3QgaGF2ZSBhdCBsZWFzdCBkKzEgcG9pbnRzXCIpXG4gIH1cbiAgdmFyIGQgPSBwb2ludHNbMF0ubGVuZ3RoXG4gIGlmKG4gPD0gZCkge1xuICAgIHRocm93IG5ldyBFcnJvcihcIk11c3QgaW5wdXQgYXQgbGVhc3QgZCsxIHBvaW50c1wiKVxuICB9XG5cbiAgLy9GSVhNRTogVGhpcyBjb3VsZCBiZSBkZWdlbmVyYXRlLCBidXQgbmVlZCB0byBzZWxlY3QgZCsxIG5vbi1jb3BsYW5hciBwb2ludHMgdG8gYm9vdHN0cmFwIHByb2Nlc3NcbiAgdmFyIGluaXRpYWxTaW1wbGV4ID0gcG9pbnRzLnNsaWNlKDAsIGQrMSlcblxuICAvL01ha2Ugc3VyZSBpbml0aWFsIHNpbXBsZXggaXMgcG9zaXRpdmVseSBvcmllbnRlZFxuICB2YXIgbyA9IG9yaWVudC5hcHBseSh2b2lkIDAsIGluaXRpYWxTaW1wbGV4KVxuICBpZihvID09PSAwKSB7XG4gICAgdGhyb3cgbmV3IEVycm9yKFwiSW5wdXQgbm90IGluIGdlbmVyYWwgcG9zaXRpb25cIilcbiAgfVxuICB2YXIgaW5pdGlhbENvb3JkcyA9IG5ldyBBcnJheShkKzEpXG4gIGZvcih2YXIgaT0wOyBpPD1kOyArK2kpIHtcbiAgICBpbml0aWFsQ29vcmRzW2ldID0gaVxuICB9XG4gIGlmKG8gPCAwKSB7XG4gICAgaW5pdGlhbENvb3Jkc1swXSA9IDFcbiAgICBpbml0aWFsQ29vcmRzWzFdID0gMFxuICB9XG5cbiAgLy9DcmVhdGUgaW5pdGlhbCB0b3BvbG9naWNhbCBpbmRleCwgZ2x1ZSBwb2ludGVycyB0b2dldGhlciAoa2luZCBvZiBtZXNzeSlcbiAgdmFyIGluaXRpYWxDZWxsID0gbmV3IFNpbXBsZXgoaW5pdGlhbENvb3JkcywgbmV3IEFycmF5KGQrMSksIGZhbHNlKVxuICB2YXIgYm91bmRhcnkgPSBpbml0aWFsQ2VsbC5hZGphY2VudFxuICB2YXIgbGlzdCA9IG5ldyBBcnJheShkKzIpXG4gIGZvcih2YXIgaT0wOyBpPD1kOyArK2kpIHtcbiAgICB2YXIgdmVydHMgPSBpbml0aWFsQ29vcmRzLnNsaWNlKClcbiAgICBmb3IodmFyIGo9MDsgajw9ZDsgKytqKSB7XG4gICAgICBpZihqID09PSBpKSB7XG4gICAgICAgIHZlcnRzW2pdID0gLTFcbiAgICAgIH1cbiAgICB9XG4gICAgdmFyIHQgPSB2ZXJ0c1swXVxuICAgIHZlcnRzWzBdID0gdmVydHNbMV1cbiAgICB2ZXJ0c1sxXSA9IHRcbiAgICB2YXIgY2VsbCA9IG5ldyBTaW1wbGV4KHZlcnRzLCBuZXcgQXJyYXkoZCsxKSwgdHJ1ZSlcbiAgICBib3VuZGFyeVtpXSA9IGNlbGxcbiAgICBsaXN0W2ldID0gY2VsbFxuICB9XG4gIGxpc3RbZCsxXSA9IGluaXRpYWxDZWxsXG4gIGZvcih2YXIgaT0wOyBpPD1kOyArK2kpIHtcbiAgICB2YXIgdmVydHMgPSBib3VuZGFyeVtpXS52ZXJ0aWNlc1xuICAgIHZhciBhZGogPSBib3VuZGFyeVtpXS5hZGphY2VudFxuICAgIGZvcih2YXIgaj0wOyBqPD1kOyArK2opIHtcbiAgICAgIHZhciB2ID0gdmVydHNbal1cbiAgICAgIGlmKHYgPCAwKSB7XG4gICAgICAgIGFkaltqXSA9IGluaXRpYWxDZWxsXG4gICAgICAgIGNvbnRpbnVlXG4gICAgICB9XG4gICAgICBmb3IodmFyIGs9MDsgazw9ZDsgKytrKSB7XG4gICAgICAgIGlmKGJvdW5kYXJ5W2tdLnZlcnRpY2VzLmluZGV4T2YodikgPCAwKSB7XG4gICAgICAgICAgYWRqW2pdID0gYm91bmRhcnlba11cbiAgICAgICAgfVxuICAgICAgfVxuICAgIH1cbiAgfVxuXG4gIC8vSW5pdGlhbGl6ZSB0cmlhbmdsZXNcbiAgdmFyIHRyaWFuZ2xlcyA9IG5ldyBUcmlhbmd1bGF0aW9uKGQsIGluaXRpYWxTaW1wbGV4LCBsaXN0KVxuXG4gIC8vSW5zZXJ0IHJlbWFpbmluZyBwb2ludHNcbiAgdmFyIHVzZVJhbmRvbSA9ICEhcmFuZG9tU2VhcmNoXG4gIGZvcih2YXIgaT1kKzE7IGk8bjsgKytpKSB7XG4gICAgdHJpYW5nbGVzLmluc2VydChwb2ludHNbaV0sIHVzZVJhbmRvbSlcbiAgfVxuICBcbiAgLy9FeHRyYWN0IGJvdW5kYXJ5IGNlbGxzXG4gIHJldHVybiB0cmlhbmdsZXMuYm91bmRhcnkoKVxufSIsIlwidXNlIHN0cmljdFwiXG5cbm1vZHVsZS5leHBvcnRzID0gZmFzdFR3b1N1bVxuXG5mdW5jdGlvbiBmYXN0VHdvU3VtKGEsIGIsIHJlc3VsdCkge1xuXHR2YXIgeCA9IGEgKyBiXG5cdHZhciBidiA9IHggLSBhXG5cdHZhciBhdiA9IHggLSBidlxuXHR2YXIgYnIgPSBiIC0gYnZcblx0dmFyIGFyID0gYSAtIGF2XG5cdGlmKHJlc3VsdCkge1xuXHRcdHJlc3VsdFswXSA9IGFyICsgYnJcblx0XHRyZXN1bHRbMV0gPSB4XG5cdFx0cmV0dXJuIHJlc3VsdFxuXHR9XG5cdHJldHVybiBbYXIrYnIsIHhdXG59IiwiXCJ1c2Ugc3RyaWN0XCJcblxudmFyIHR3b1Byb2R1Y3QgPSByZXF1aXJlKFwidHdvLXByb2R1Y3RcIilcbnZhciB0d29TdW0gPSByZXF1aXJlKFwidHdvLXN1bVwiKVxuXG5tb2R1bGUuZXhwb3J0cyA9IHNjYWxlTGluZWFyRXhwYW5zaW9uXG5cbmZ1bmN0aW9uIHNjYWxlTGluZWFyRXhwYW5zaW9uKGUsIHNjYWxlKSB7XG4gIHZhciBuID0gZS5sZW5ndGhcbiAgaWYobiA9PT0gMSkge1xuICAgIHZhciB0cyA9IHR3b1Byb2R1Y3QoZVswXSwgc2NhbGUpXG4gICAgaWYodHNbMF0pIHtcbiAgICAgIHJldHVybiB0c1xuICAgIH1cbiAgICByZXR1cm4gWyB0c1sxXSBdXG4gIH1cbiAgdmFyIGcgPSBuZXcgQXJyYXkoMiAqIG4pXG4gIHZhciBxID0gWzAuMSwgMC4xXVxuICB2YXIgdCA9IFswLjEsIDAuMV1cbiAgdmFyIGNvdW50ID0gMFxuICB0d29Qcm9kdWN0KGVbMF0sIHNjYWxlLCBxKVxuICBpZihxWzBdKSB7XG4gICAgZ1tjb3VudCsrXSA9IHFbMF1cbiAgfVxuICBmb3IodmFyIGk9MTsgaTxuOyArK2kpIHtcbiAgICB0d29Qcm9kdWN0KGVbaV0sIHNjYWxlLCB0KVxuICAgIHZhciBwcSA9IHFbMV1cbiAgICB0d29TdW0ocHEsIHRbMF0sIHEpXG4gICAgaWYocVswXSkge1xuICAgICAgZ1tjb3VudCsrXSA9IHFbMF1cbiAgICB9XG4gICAgdmFyIGEgPSB0WzFdXG4gICAgdmFyIGIgPSBxWzFdXG4gICAgdmFyIHggPSBhICsgYlxuICAgIHZhciBidiA9IHggLSBhXG4gICAgdmFyIHkgPSBiIC0gYnZcbiAgICBxWzFdID0geFxuICAgIGlmKHkpIHtcbiAgICAgIGdbY291bnQrK10gPSB5XG4gICAgfVxuICB9XG4gIGlmKHFbMV0pIHtcbiAgICBnW2NvdW50KytdID0gcVsxXVxuICB9XG4gIGlmKGNvdW50ID09PSAwKSB7XG4gICAgZ1tjb3VudCsrXSA9IDAuMFxuICB9XG4gIGcubGVuZ3RoID0gY291bnRcbiAgcmV0dXJuIGdcbn0iLCJcInVzZSBzdHJpY3RcIlxuXG5tb2R1bGUuZXhwb3J0cyA9IHJvYnVzdFN1YnRyYWN0XG5cbi8vRWFzeSBjYXNlOiBBZGQgdHdvIHNjYWxhcnNcbmZ1bmN0aW9uIHNjYWxhclNjYWxhcihhLCBiKSB7XG4gIHZhciB4ID0gYSArIGJcbiAgdmFyIGJ2ID0geCAtIGFcbiAgdmFyIGF2ID0geCAtIGJ2XG4gIHZhciBiciA9IGIgLSBidlxuICB2YXIgYXIgPSBhIC0gYXZcbiAgdmFyIHkgPSBhciArIGJyXG4gIGlmKHkpIHtcbiAgICByZXR1cm4gW3ksIHhdXG4gIH1cbiAgcmV0dXJuIFt4XVxufVxuXG5mdW5jdGlvbiByb2J1c3RTdWJ0cmFjdChlLCBmKSB7XG4gIHZhciBuZSA9IGUubGVuZ3RofDBcbiAgdmFyIG5mID0gZi5sZW5ndGh8MFxuICBpZihuZSA9PT0gMSAmJiBuZiA9PT0gMSkge1xuICAgIHJldHVybiBzY2FsYXJTY2FsYXIoZVswXSwgLWZbMF0pXG4gIH1cbiAgdmFyIG4gPSBuZSArIG5mXG4gIHZhciBnID0gbmV3IEFycmF5KG4pXG4gIHZhciBjb3VudCA9IDBcbiAgdmFyIGVwdHIgPSAwXG4gIHZhciBmcHRyID0gMFxuICB2YXIgYWJzID0gTWF0aC5hYnNcbiAgdmFyIGVpID0gZVtlcHRyXVxuICB2YXIgZWEgPSBhYnMoZWkpXG4gIHZhciBmaSA9IC1mW2ZwdHJdXG4gIHZhciBmYSA9IGFicyhmaSlcbiAgdmFyIGEsIGJcbiAgaWYoZWEgPCBmYSkge1xuICAgIGIgPSBlaVxuICAgIGVwdHIgKz0gMVxuICAgIGlmKGVwdHIgPCBuZSkge1xuICAgICAgZWkgPSBlW2VwdHJdXG4gICAgICBlYSA9IGFicyhlaSlcbiAgICB9XG4gIH0gZWxzZSB7XG4gICAgYiA9IGZpXG4gICAgZnB0ciArPSAxXG4gICAgaWYoZnB0ciA8IG5mKSB7XG4gICAgICBmaSA9IC1mW2ZwdHJdXG4gICAgICBmYSA9IGFicyhmaSlcbiAgICB9XG4gIH1cbiAgaWYoKGVwdHIgPCBuZSAmJiBlYSA8IGZhKSB8fCAoZnB0ciA+PSBuZikpIHtcbiAgICBhID0gZWlcbiAgICBlcHRyICs9IDFcbiAgICBpZihlcHRyIDwgbmUpIHtcbiAgICAgIGVpID0gZVtlcHRyXVxuICAgICAgZWEgPSBhYnMoZWkpXG4gICAgfVxuICB9IGVsc2Uge1xuICAgIGEgPSBmaVxuICAgIGZwdHIgKz0gMVxuICAgIGlmKGZwdHIgPCBuZikge1xuICAgICAgZmkgPSAtZltmcHRyXVxuICAgICAgZmEgPSBhYnMoZmkpXG4gICAgfVxuICB9XG4gIHZhciB4ID0gYSArIGJcbiAgdmFyIGJ2ID0geCAtIGFcbiAgdmFyIHkgPSBiIC0gYnZcbiAgdmFyIHEwID0geVxuICB2YXIgcTEgPSB4XG4gIHZhciBfeCwgX2J2LCBfYXYsIF9iciwgX2FyXG4gIHdoaWxlKGVwdHIgPCBuZSAmJiBmcHRyIDwgbmYpIHtcbiAgICBpZihlYSA8IGZhKSB7XG4gICAgICBhID0gZWlcbiAgICAgIGVwdHIgKz0gMVxuICAgICAgaWYoZXB0ciA8IG5lKSB7XG4gICAgICAgIGVpID0gZVtlcHRyXVxuICAgICAgICBlYSA9IGFicyhlaSlcbiAgICAgIH1cbiAgICB9IGVsc2Uge1xuICAgICAgYSA9IGZpXG4gICAgICBmcHRyICs9IDFcbiAgICAgIGlmKGZwdHIgPCBuZikge1xuICAgICAgICBmaSA9IC1mW2ZwdHJdXG4gICAgICAgIGZhID0gYWJzKGZpKVxuICAgICAgfVxuICAgIH1cbiAgICBiID0gcTBcbiAgICB4ID0gYSArIGJcbiAgICBidiA9IHggLSBhXG4gICAgeSA9IGIgLSBidlxuICAgIGlmKHkpIHtcbiAgICAgIGdbY291bnQrK10gPSB5XG4gICAgfVxuICAgIF94ID0gcTEgKyB4XG4gICAgX2J2ID0gX3ggLSBxMVxuICAgIF9hdiA9IF94IC0gX2J2XG4gICAgX2JyID0geCAtIF9idlxuICAgIF9hciA9IHExIC0gX2F2XG4gICAgcTAgPSBfYXIgKyBfYnJcbiAgICBxMSA9IF94XG4gIH1cbiAgd2hpbGUoZXB0ciA8IG5lKSB7XG4gICAgYSA9IGVpXG4gICAgYiA9IHEwXG4gICAgeCA9IGEgKyBiXG4gICAgYnYgPSB4IC0gYVxuICAgIHkgPSBiIC0gYnZcbiAgICBpZih5KSB7XG4gICAgICBnW2NvdW50KytdID0geVxuICAgIH1cbiAgICBfeCA9IHExICsgeFxuICAgIF9idiA9IF94IC0gcTFcbiAgICBfYXYgPSBfeCAtIF9idlxuICAgIF9iciA9IHggLSBfYnZcbiAgICBfYXIgPSBxMSAtIF9hdlxuICAgIHEwID0gX2FyICsgX2JyXG4gICAgcTEgPSBfeFxuICAgIGVwdHIgKz0gMVxuICAgIGlmKGVwdHIgPCBuZSkge1xuICAgICAgZWkgPSBlW2VwdHJdXG4gICAgfVxuICB9XG4gIHdoaWxlKGZwdHIgPCBuZikge1xuICAgIGEgPSBmaVxuICAgIGIgPSBxMFxuICAgIHggPSBhICsgYlxuICAgIGJ2ID0geCAtIGFcbiAgICB5ID0gYiAtIGJ2XG4gICAgaWYoeSkge1xuICAgICAgZ1tjb3VudCsrXSA9IHlcbiAgICB9IFxuICAgIF94ID0gcTEgKyB4XG4gICAgX2J2ID0gX3ggLSBxMVxuICAgIF9hdiA9IF94IC0gX2J2XG4gICAgX2JyID0geCAtIF9idlxuICAgIF9hciA9IHExIC0gX2F2XG4gICAgcTAgPSBfYXIgKyBfYnJcbiAgICBxMSA9IF94XG4gICAgZnB0ciArPSAxXG4gICAgaWYoZnB0ciA8IG5mKSB7XG4gICAgICBmaSA9IC1mW2ZwdHJdXG4gICAgfVxuICB9XG4gIGlmKHEwKSB7XG4gICAgZ1tjb3VudCsrXSA9IHEwXG4gIH1cbiAgaWYocTEpIHtcbiAgICBnW2NvdW50KytdID0gcTFcbiAgfVxuICBpZighY291bnQpIHtcbiAgICBnW2NvdW50KytdID0gMC4wICBcbiAgfVxuICBnLmxlbmd0aCA9IGNvdW50XG4gIHJldHVybiBnXG59IiwiXCJ1c2Ugc3RyaWN0XCJcblxubW9kdWxlLmV4cG9ydHMgPSBsaW5lYXJFeHBhbnNpb25TdW1cblxuLy9FYXN5IGNhc2U6IEFkZCB0d28gc2NhbGFyc1xuZnVuY3Rpb24gc2NhbGFyU2NhbGFyKGEsIGIpIHtcbiAgdmFyIHggPSBhICsgYlxuICB2YXIgYnYgPSB4IC0gYVxuICB2YXIgYXYgPSB4IC0gYnZcbiAgdmFyIGJyID0gYiAtIGJ2XG4gIHZhciBhciA9IGEgLSBhdlxuICB2YXIgeSA9IGFyICsgYnJcbiAgaWYoeSkge1xuICAgIHJldHVybiBbeSwgeF1cbiAgfVxuICByZXR1cm4gW3hdXG59XG5cbmZ1bmN0aW9uIGxpbmVhckV4cGFuc2lvblN1bShlLCBmKSB7XG4gIHZhciBuZSA9IGUubGVuZ3RofDBcbiAgdmFyIG5mID0gZi5sZW5ndGh8MFxuICBpZihuZSA9PT0gMSAmJiBuZiA9PT0gMSkge1xuICAgIHJldHVybiBzY2FsYXJTY2FsYXIoZVswXSwgZlswXSlcbiAgfVxuICB2YXIgbiA9IG5lICsgbmZcbiAgdmFyIGcgPSBuZXcgQXJyYXkobilcbiAgdmFyIGNvdW50ID0gMFxuICB2YXIgZXB0ciA9IDBcbiAgdmFyIGZwdHIgPSAwXG4gIHZhciBhYnMgPSBNYXRoLmFic1xuICB2YXIgZWkgPSBlW2VwdHJdXG4gIHZhciBlYSA9IGFicyhlaSlcbiAgdmFyIGZpID0gZltmcHRyXVxuICB2YXIgZmEgPSBhYnMoZmkpXG4gIHZhciBhLCBiXG4gIGlmKGVhIDwgZmEpIHtcbiAgICBiID0gZWlcbiAgICBlcHRyICs9IDFcbiAgICBpZihlcHRyIDwgbmUpIHtcbiAgICAgIGVpID0gZVtlcHRyXVxuICAgICAgZWEgPSBhYnMoZWkpXG4gICAgfVxuICB9IGVsc2Uge1xuICAgIGIgPSBmaVxuICAgIGZwdHIgKz0gMVxuICAgIGlmKGZwdHIgPCBuZikge1xuICAgICAgZmkgPSBmW2ZwdHJdXG4gICAgICBmYSA9IGFicyhmaSlcbiAgICB9XG4gIH1cbiAgaWYoKGVwdHIgPCBuZSAmJiBlYSA8IGZhKSB8fCAoZnB0ciA+PSBuZikpIHtcbiAgICBhID0gZWlcbiAgICBlcHRyICs9IDFcbiAgICBpZihlcHRyIDwgbmUpIHtcbiAgICAgIGVpID0gZVtlcHRyXVxuICAgICAgZWEgPSBhYnMoZWkpXG4gICAgfVxuICB9IGVsc2Uge1xuICAgIGEgPSBmaVxuICAgIGZwdHIgKz0gMVxuICAgIGlmKGZwdHIgPCBuZikge1xuICAgICAgZmkgPSBmW2ZwdHJdXG4gICAgICBmYSA9IGFicyhmaSlcbiAgICB9XG4gIH1cbiAgdmFyIHggPSBhICsgYlxuICB2YXIgYnYgPSB4IC0gYVxuICB2YXIgeSA9IGIgLSBidlxuICB2YXIgcTAgPSB5XG4gIHZhciBxMSA9IHhcbiAgdmFyIF94LCBfYnYsIF9hdiwgX2JyLCBfYXJcbiAgd2hpbGUoZXB0ciA8IG5lICYmIGZwdHIgPCBuZikge1xuICAgIGlmKGVhIDwgZmEpIHtcbiAgICAgIGEgPSBlaVxuICAgICAgZXB0ciArPSAxXG4gICAgICBpZihlcHRyIDwgbmUpIHtcbiAgICAgICAgZWkgPSBlW2VwdHJdXG4gICAgICAgIGVhID0gYWJzKGVpKVxuICAgICAgfVxuICAgIH0gZWxzZSB7XG4gICAgICBhID0gZmlcbiAgICAgIGZwdHIgKz0gMVxuICAgICAgaWYoZnB0ciA8IG5mKSB7XG4gICAgICAgIGZpID0gZltmcHRyXVxuICAgICAgICBmYSA9IGFicyhmaSlcbiAgICAgIH1cbiAgICB9XG4gICAgYiA9IHEwXG4gICAgeCA9IGEgKyBiXG4gICAgYnYgPSB4IC0gYVxuICAgIHkgPSBiIC0gYnZcbiAgICBpZih5KSB7XG4gICAgICBnW2NvdW50KytdID0geVxuICAgIH1cbiAgICBfeCA9IHExICsgeFxuICAgIF9idiA9IF94IC0gcTFcbiAgICBfYXYgPSBfeCAtIF9idlxuICAgIF9iciA9IHggLSBfYnZcbiAgICBfYXIgPSBxMSAtIF9hdlxuICAgIHEwID0gX2FyICsgX2JyXG4gICAgcTEgPSBfeFxuICB9XG4gIHdoaWxlKGVwdHIgPCBuZSkge1xuICAgIGEgPSBlaVxuICAgIGIgPSBxMFxuICAgIHggPSBhICsgYlxuICAgIGJ2ID0geCAtIGFcbiAgICB5ID0gYiAtIGJ2XG4gICAgaWYoeSkge1xuICAgICAgZ1tjb3VudCsrXSA9IHlcbiAgICB9XG4gICAgX3ggPSBxMSArIHhcbiAgICBfYnYgPSBfeCAtIHExXG4gICAgX2F2ID0gX3ggLSBfYnZcbiAgICBfYnIgPSB4IC0gX2J2XG4gICAgX2FyID0gcTEgLSBfYXZcbiAgICBxMCA9IF9hciArIF9iclxuICAgIHExID0gX3hcbiAgICBlcHRyICs9IDFcbiAgICBpZihlcHRyIDwgbmUpIHtcbiAgICAgIGVpID0gZVtlcHRyXVxuICAgIH1cbiAgfVxuICB3aGlsZShmcHRyIDwgbmYpIHtcbiAgICBhID0gZmlcbiAgICBiID0gcTBcbiAgICB4ID0gYSArIGJcbiAgICBidiA9IHggLSBhXG4gICAgeSA9IGIgLSBidlxuICAgIGlmKHkpIHtcbiAgICAgIGdbY291bnQrK10gPSB5XG4gICAgfSBcbiAgICBfeCA9IHExICsgeFxuICAgIF9idiA9IF94IC0gcTFcbiAgICBfYXYgPSBfeCAtIF9idlxuICAgIF9iciA9IHggLSBfYnZcbiAgICBfYXIgPSBxMSAtIF9hdlxuICAgIHEwID0gX2FyICsgX2JyXG4gICAgcTEgPSBfeFxuICAgIGZwdHIgKz0gMVxuICAgIGlmKGZwdHIgPCBuZikge1xuICAgICAgZmkgPSBmW2ZwdHJdXG4gICAgfVxuICB9XG4gIGlmKHEwKSB7XG4gICAgZ1tjb3VudCsrXSA9IHEwXG4gIH1cbiAgaWYocTEpIHtcbiAgICBnW2NvdW50KytdID0gcTFcbiAgfVxuICBpZighY291bnQpIHtcbiAgICBnW2NvdW50KytdID0gMC4wICBcbiAgfVxuICBnLmxlbmd0aCA9IGNvdW50XG4gIHJldHVybiBnXG59IiwiXCJ1c2Ugc3RyaWN0XCJcblxubW9kdWxlLmV4cG9ydHMgPSB0d29Qcm9kdWN0XG5cbnZhciBTUExJVFRFUiA9ICsoTWF0aC5wb3coMiwgMjcpICsgMS4wKVxuXG5mdW5jdGlvbiB0d29Qcm9kdWN0KGEsIGIsIHJlc3VsdCkge1xuICB2YXIgeCA9IGEgKiBiXG5cbiAgdmFyIGMgPSBTUExJVFRFUiAqIGFcbiAgdmFyIGFiaWcgPSBjIC0gYVxuICB2YXIgYWhpID0gYyAtIGFiaWdcbiAgdmFyIGFsbyA9IGEgLSBhaGlcblxuICB2YXIgZCA9IFNQTElUVEVSICogYlxuICB2YXIgYmJpZyA9IGQgLSBiXG4gIHZhciBiaGkgPSBkIC0gYmJpZ1xuICB2YXIgYmxvID0gYiAtIGJoaVxuXG4gIHZhciBlcnIxID0geCAtIChhaGkgKiBiaGkpXG4gIHZhciBlcnIyID0gZXJyMSAtIChhbG8gKiBiaGkpXG4gIHZhciBlcnIzID0gZXJyMiAtIChhaGkgKiBibG8pXG5cbiAgdmFyIHkgPSBhbG8gKiBibG8gLSBlcnIzXG5cbiAgaWYocmVzdWx0KSB7XG4gICAgcmVzdWx0WzBdID0geVxuICAgIHJlc3VsdFsxXSA9IHhcbiAgICByZXR1cm4gcmVzdWx0XG4gIH1cblxuICByZXR1cm4gWyB5LCB4IF1cbn0iLCJcInVzZSBzdHJpY3RcIlxuXG52YXIgdHdvUHJvZHVjdCA9IHJlcXVpcmUoXCJ0d28tcHJvZHVjdFwiKVxudmFyIHJvYnVzdFN1bSA9IHJlcXVpcmUoXCJyb2J1c3Qtc3VtXCIpXG52YXIgcm9idXN0U2NhbGUgPSByZXF1aXJlKFwicm9idXN0LXNjYWxlXCIpXG52YXIgcm9idXN0U3VidHJhY3QgPSByZXF1aXJlKFwicm9idXN0LXN1YnRyYWN0XCIpXG5cbnZhciBOVU1fRVhQQU5EID0gNVxuXG52YXIgRVBTSUxPTiAgICAgPSAxLjExMDIyMzAyNDYyNTE1NjVlLTE2XG52YXIgRVJSQk9VTkQzICAgPSAoMy4wICsgMTYuMCAqIEVQU0lMT04pICogRVBTSUxPTlxudmFyIEVSUkJPVU5ENCAgID0gKDcuMCArIDU2LjAgKiBFUFNJTE9OKSAqIEVQU0lMT05cblxuZnVuY3Rpb24gY29mYWN0b3IobSwgYykge1xuICB2YXIgcmVzdWx0ID0gbmV3IEFycmF5KG0ubGVuZ3RoLTEpXG4gIGZvcih2YXIgaT0xOyBpPG0ubGVuZ3RoOyArK2kpIHtcbiAgICB2YXIgciA9IHJlc3VsdFtpLTFdID0gbmV3IEFycmF5KG0ubGVuZ3RoLTEpXG4gICAgZm9yKHZhciBqPTAsaz0wOyBqPG0ubGVuZ3RoOyArK2opIHtcbiAgICAgIGlmKGogPT09IGMpIHtcbiAgICAgICAgY29udGludWVcbiAgICAgIH1cbiAgICAgIHJbaysrXSA9IG1baV1bal1cbiAgICB9XG4gIH1cbiAgcmV0dXJuIHJlc3VsdFxufVxuXG5mdW5jdGlvbiBtYXRyaXgobikge1xuICB2YXIgcmVzdWx0ID0gbmV3IEFycmF5KG4pXG4gIGZvcih2YXIgaT0wOyBpPG47ICsraSkge1xuICAgIHJlc3VsdFtpXSA9IG5ldyBBcnJheShuKVxuICAgIGZvcih2YXIgaj0wOyBqPG47ICsraikge1xuICAgICAgcmVzdWx0W2ldW2pdID0gW1wibVwiLCBqLCBcIltcIiwgKG4taS0xKSwgXCJdXCJdLmpvaW4oXCJcIilcbiAgICB9XG4gIH1cbiAgcmV0dXJuIHJlc3VsdFxufVxuXG5mdW5jdGlvbiBzaWduKG4pIHtcbiAgaWYobiAmIDEpIHtcbiAgICByZXR1cm4gXCItXCJcbiAgfVxuICByZXR1cm4gXCJcIlxufVxuXG5mdW5jdGlvbiBnZW5lcmF0ZVN1bShleHByKSB7XG4gIGlmKGV4cHIubGVuZ3RoID09PSAxKSB7XG4gICAgcmV0dXJuIGV4cHJbMF1cbiAgfSBlbHNlIGlmKGV4cHIubGVuZ3RoID09PSAyKSB7XG4gICAgcmV0dXJuIFtcInN1bShcIiwgZXhwclswXSwgXCIsXCIsIGV4cHJbMV0sIFwiKVwiXS5qb2luKFwiXCIpXG4gIH0gZWxzZSB7XG4gICAgdmFyIG0gPSBleHByLmxlbmd0aD4+MVxuICAgIHJldHVybiBbXCJzdW0oXCIsIGdlbmVyYXRlU3VtKGV4cHIuc2xpY2UoMCwgbSkpLCBcIixcIiwgZ2VuZXJhdGVTdW0oZXhwci5zbGljZShtKSksIFwiKVwiXS5qb2luKFwiXCIpXG4gIH1cbn1cblxuZnVuY3Rpb24gZGV0ZXJtaW5hbnQobSkge1xuICBpZihtLmxlbmd0aCA9PT0gMikge1xuICAgIHJldHVybiBbW1wic3VtKHByb2QoXCIsIG1bMF1bMF0sIFwiLFwiLCBtWzFdWzFdLCBcIikscHJvZCgtXCIsIG1bMF1bMV0sIFwiLFwiLCBtWzFdWzBdLCBcIikpXCJdLmpvaW4oXCJcIildXG4gIH0gZWxzZSB7XG4gICAgdmFyIGV4cHIgPSBbXVxuICAgIGZvcih2YXIgaT0wOyBpPG0ubGVuZ3RoOyArK2kpIHtcbiAgICAgIGV4cHIucHVzaChbXCJzY2FsZShcIiwgZ2VuZXJhdGVTdW0oZGV0ZXJtaW5hbnQoY29mYWN0b3IobSwgaSkpKSwgXCIsXCIsIHNpZ24oaSksIG1bMF1baV0sIFwiKVwiXS5qb2luKFwiXCIpKVxuICAgIH1cbiAgICByZXR1cm4gZXhwclxuICB9XG59XG5cbmZ1bmN0aW9uIG9yaWVudGF0aW9uKG4pIHtcbiAgdmFyIHBvcyA9IFtdXG4gIHZhciBuZWcgPSBbXVxuICB2YXIgbSA9IG1hdHJpeChuKVxuICB2YXIgYXJncyA9IFtdXG4gIGZvcih2YXIgaT0wOyBpPG47ICsraSkge1xuICAgIGlmKChpJjEpPT09MCkge1xuICAgICAgcG9zLnB1c2guYXBwbHkocG9zLCBkZXRlcm1pbmFudChjb2ZhY3RvcihtLCBpKSkpXG4gICAgfSBlbHNlIHtcbiAgICAgIG5lZy5wdXNoLmFwcGx5KG5lZywgZGV0ZXJtaW5hbnQoY29mYWN0b3IobSwgaSkpKVxuICAgIH1cbiAgICBhcmdzLnB1c2goXCJtXCIgKyBpKVxuICB9XG4gIHZhciBwb3NFeHByID0gZ2VuZXJhdGVTdW0ocG9zKVxuICB2YXIgbmVnRXhwciA9IGdlbmVyYXRlU3VtKG5lZylcbiAgdmFyIGZ1bmNOYW1lID0gXCJvcmllbnRhdGlvblwiICsgbiArIFwiRXhhY3RcIlxuICB2YXIgY29kZSA9IFtcImZ1bmN0aW9uIFwiLCBmdW5jTmFtZSwgXCIoXCIsIGFyZ3Muam9pbigpLCBcIil7dmFyIHA9XCIsIHBvc0V4cHIsIFwiLG49XCIsIG5lZ0V4cHIsIFwiLGQ9c3ViKHAsbik7XFxcbnJldHVybiBkW2QubGVuZ3RoLTFdO307cmV0dXJuIFwiLCBmdW5jTmFtZV0uam9pbihcIlwiKVxuICB2YXIgcHJvYyA9IG5ldyBGdW5jdGlvbihcInN1bVwiLCBcInByb2RcIiwgXCJzY2FsZVwiLCBcInN1YlwiLCBjb2RlKVxuICByZXR1cm4gcHJvYyhyb2J1c3RTdW0sIHR3b1Byb2R1Y3QsIHJvYnVzdFNjYWxlLCByb2J1c3RTdWJ0cmFjdClcbn1cblxudmFyIG9yaWVudGF0aW9uM0V4YWN0ID0gb3JpZW50YXRpb24oMylcbnZhciBvcmllbnRhdGlvbjRFeGFjdCA9IG9yaWVudGF0aW9uKDQpXG5cbnZhciBDQUNIRUQgPSBbXG4gIGZ1bmN0aW9uIG9yaWVudGF0aW9uMCgpIHsgcmV0dXJuIDAgfSxcbiAgZnVuY3Rpb24gb3JpZW50YXRpb24xKCkgeyByZXR1cm4gMCB9LFxuICBmdW5jdGlvbiBvcmllbnRhdGlvbjIoYSwgYikgeyBcbiAgICByZXR1cm4gYlswXSAtIGFbMF1cbiAgfSxcbiAgZnVuY3Rpb24gb3JpZW50YXRpb24zKGEsIGIsIGMpIHtcbiAgICB2YXIgbCA9IChhWzFdIC0gY1sxXSkgKiAoYlswXSAtIGNbMF0pXG4gICAgdmFyIHIgPSAoYVswXSAtIGNbMF0pICogKGJbMV0gLSBjWzFdKVxuICAgIHZhciBkZXQgPSBsIC0gclxuICAgIHZhciBzXG4gICAgaWYobCA+IDApIHtcbiAgICAgIGlmKHIgPD0gMCkge1xuICAgICAgICByZXR1cm4gZGV0XG4gICAgICB9IGVsc2Uge1xuICAgICAgICBzID0gbCArIHJcbiAgICAgIH1cbiAgICB9IGVsc2UgaWYobCA8IDApIHtcbiAgICAgIGlmKHIgPj0gMCkge1xuICAgICAgICByZXR1cm4gZGV0XG4gICAgICB9IGVsc2Uge1xuICAgICAgICBzID0gLShsICsgcilcbiAgICAgIH1cbiAgICB9IGVsc2Uge1xuICAgICAgcmV0dXJuIGRldFxuICAgIH1cbiAgICB2YXIgdG9sID0gRVJSQk9VTkQzICogc1xuICAgIGlmKGRldCA+PSB0b2wgfHwgZGV0IDw9IC10b2wpIHtcbiAgICAgIHJldHVybiBkZXRcbiAgICB9XG4gICAgcmV0dXJuIG9yaWVudGF0aW9uM0V4YWN0KGEsIGIsIGMpXG4gIH0sXG4gIGZ1bmN0aW9uIG9yaWVudGF0aW9uNChhLGIsYyxkKSB7XG4gICAgdmFyIGFkeCA9IGFbMF0gLSBkWzBdXG4gICAgdmFyIGJkeCA9IGJbMF0gLSBkWzBdXG4gICAgdmFyIGNkeCA9IGNbMF0gLSBkWzBdXG4gICAgdmFyIGFkeSA9IGFbMV0gLSBkWzFdXG4gICAgdmFyIGJkeSA9IGJbMV0gLSBkWzFdXG4gICAgdmFyIGNkeSA9IGNbMV0gLSBkWzFdXG4gICAgdmFyIGFkeiA9IGFbMl0gLSBkWzJdXG4gICAgdmFyIGJkeiA9IGJbMl0gLSBkWzJdXG4gICAgdmFyIGNkeiA9IGNbMl0gLSBkWzJdXG4gICAgdmFyIGJkeGNkeSA9IGJkeCAqIGNkeVxuICAgIHZhciBjZHhiZHkgPSBjZHggKiBiZHlcbiAgICB2YXIgY2R4YWR5ID0gY2R4ICogYWR5XG4gICAgdmFyIGFkeGNkeSA9IGFkeCAqIGNkeVxuICAgIHZhciBhZHhiZHkgPSBhZHggKiBiZHlcbiAgICB2YXIgYmR4YWR5ID0gYmR4ICogYWR5XG4gICAgdmFyIGRldCA9IGFkeiAqIChiZHhjZHkgLSBjZHhiZHkpIFxuICAgICAgICAgICAgKyBiZHogKiAoY2R4YWR5IC0gYWR4Y2R5KVxuICAgICAgICAgICAgKyBjZHogKiAoYWR4YmR5IC0gYmR4YWR5KVxuICAgIHZhciBwZXJtYW5lbnQgPSAoTWF0aC5hYnMoYmR4Y2R5KSArIE1hdGguYWJzKGNkeGJkeSkpICogTWF0aC5hYnMoYWR6KVxuICAgICAgICAgICAgICAgICAgKyAoTWF0aC5hYnMoY2R4YWR5KSArIE1hdGguYWJzKGFkeGNkeSkpICogTWF0aC5hYnMoYmR6KVxuICAgICAgICAgICAgICAgICAgKyAoTWF0aC5hYnMoYWR4YmR5KSArIE1hdGguYWJzKGJkeGFkeSkpICogTWF0aC5hYnMoY2R6KVxuICAgIHZhciB0b2wgPSBFUlJCT1VORDQgKiBwZXJtYW5lbnRcbiAgICBpZiAoKGRldCA+IHRvbCkgfHwgKC1kZXQgPiB0b2wpKSB7XG4gICAgICByZXR1cm4gZGV0XG4gICAgfVxuICAgIHJldHVybiBvcmllbnRhdGlvbjRFeGFjdChhLGIsYyxkKVxuICB9XG5dXG5cbmZ1bmN0aW9uIHNsb3dPcmllbnQoYXJncykge1xuICB2YXIgcHJvYyA9IENBQ0hFRFthcmdzLmxlbmd0aF1cbiAgaWYoIXByb2MpIHtcbiAgICBwcm9jID0gQ0FDSEVEW2FyZ3MubGVuZ3RoXSA9IG9yaWVudGF0aW9uKGFyZ3MubGVuZ3RoKVxuICB9XG4gIHJldHVybiBwcm9jLmFwcGx5KHVuZGVmaW5lZCwgYXJncylcbn1cblxuZnVuY3Rpb24gZ2VuZXJhdGVPcmllbnRhdGlvblByb2MoKSB7XG4gIHdoaWxlKENBQ0hFRC5sZW5ndGggPD0gTlVNX0VYUEFORCkge1xuICAgIENBQ0hFRC5wdXNoKG9yaWVudGF0aW9uKENBQ0hFRC5sZW5ndGgpKVxuICB9XG4gIHZhciBhcmdzID0gW11cbiAgdmFyIHByb2NBcmdzID0gW1wic2xvd1wiXVxuICBmb3IodmFyIGk9MDsgaTw9TlVNX0VYUEFORDsgKytpKSB7XG4gICAgYXJncy5wdXNoKFwiYVwiICsgaSlcbiAgICBwcm9jQXJncy5wdXNoKFwib1wiICsgaSlcbiAgfVxuICB2YXIgY29kZSA9IFtcbiAgICBcImZ1bmN0aW9uIGdldE9yaWVudGF0aW9uKFwiLCBhcmdzLmpvaW4oKSwgXCIpe3N3aXRjaChhcmd1bWVudHMubGVuZ3RoKXtjYXNlIDA6Y2FzZSAxOnJldHVybiAwO1wiXG4gIF1cbiAgZm9yKHZhciBpPTI7IGk8PU5VTV9FWFBBTkQ7ICsraSkge1xuICAgIGNvZGUucHVzaChcImNhc2UgXCIsIGksIFwiOnJldHVybiBvXCIsIGksIFwiKFwiLCBhcmdzLnNsaWNlKDAsIGkpLmpvaW4oKSwgXCIpO1wiKVxuICB9XG4gIGNvZGUucHVzaChcIn12YXIgcz1uZXcgQXJyYXkoYXJndW1lbnRzLmxlbmd0aCk7Zm9yKHZhciBpPTA7aTxhcmd1bWVudHMubGVuZ3RoOysraSl7c1tpXT1hcmd1bWVudHNbaV19O3JldHVybiBzbG93KHMpO31yZXR1cm4gZ2V0T3JpZW50YXRpb25cIilcbiAgcHJvY0FyZ3MucHVzaChjb2RlLmpvaW4oXCJcIikpXG5cbiAgdmFyIHByb2MgPSBGdW5jdGlvbi5hcHBseSh1bmRlZmluZWQsIHByb2NBcmdzKVxuICBtb2R1bGUuZXhwb3J0cyA9IHByb2MuYXBwbHkodW5kZWZpbmVkLCBbc2xvd09yaWVudF0uY29uY2F0KENBQ0hFRCkpXG4gIGZvcih2YXIgaT0wOyBpPD1OVU1fRVhQQU5EOyArK2kpIHtcbiAgICBtb2R1bGUuZXhwb3J0c1tpXSA9IENBQ0hFRFtpXVxuICB9XG59XG5cbmdlbmVyYXRlT3JpZW50YXRpb25Qcm9jKCkiLCIvKipcbiAqIEJpdCB0d2lkZGxpbmcgaGFja3MgZm9yIEphdmFTY3JpcHQuXG4gKlxuICogQXV0aG9yOiBNaWtvbGEgTHlzZW5rb1xuICpcbiAqIFBvcnRlZCBmcm9tIFN0YW5mb3JkIGJpdCB0d2lkZGxpbmcgaGFjayBsaWJyYXJ5OlxuICogICAgaHR0cDovL2dyYXBoaWNzLnN0YW5mb3JkLmVkdS9+c2VhbmRlci9iaXRoYWNrcy5odG1sXG4gKi9cblxuXCJ1c2Ugc3RyaWN0XCI7IFwidXNlIHJlc3RyaWN0XCI7XG5cbi8vTnVtYmVyIG9mIGJpdHMgaW4gYW4gaW50ZWdlclxudmFyIElOVF9CSVRTID0gMzI7XG5cbi8vQ29uc3RhbnRzXG5leHBvcnRzLklOVF9CSVRTICA9IElOVF9CSVRTO1xuZXhwb3J0cy5JTlRfTUFYICAgPSAgMHg3ZmZmZmZmZjtcbmV4cG9ydHMuSU5UX01JTiAgID0gLTE8PChJTlRfQklUUy0xKTtcblxuLy9SZXR1cm5zIC0xLCAwLCArMSBkZXBlbmRpbmcgb24gc2lnbiBvZiB4XG5leHBvcnRzLnNpZ24gPSBmdW5jdGlvbih2KSB7XG4gIHJldHVybiAodiA+IDApIC0gKHYgPCAwKTtcbn1cblxuLy9Db21wdXRlcyBhYnNvbHV0ZSB2YWx1ZSBvZiBpbnRlZ2VyXG5leHBvcnRzLmFicyA9IGZ1bmN0aW9uKHYpIHtcbiAgdmFyIG1hc2sgPSB2ID4+IChJTlRfQklUUy0xKTtcbiAgcmV0dXJuICh2IF4gbWFzaykgLSBtYXNrO1xufVxuXG4vL0NvbXB1dGVzIG1pbmltdW0gb2YgaW50ZWdlcnMgeCBhbmQgeVxuZXhwb3J0cy5taW4gPSBmdW5jdGlvbih4LCB5KSB7XG4gIHJldHVybiB5IF4gKCh4IF4geSkgJiAtKHggPCB5KSk7XG59XG5cbi8vQ29tcHV0ZXMgbWF4aW11bSBvZiBpbnRlZ2VycyB4IGFuZCB5XG5leHBvcnRzLm1heCA9IGZ1bmN0aW9uKHgsIHkpIHtcbiAgcmV0dXJuIHggXiAoKHggXiB5KSAmIC0oeCA8IHkpKTtcbn1cblxuLy9DaGVja3MgaWYgYSBudW1iZXIgaXMgYSBwb3dlciBvZiB0d29cbmV4cG9ydHMuaXNQb3cyID0gZnVuY3Rpb24odikge1xuICByZXR1cm4gISh2ICYgKHYtMSkpICYmICghIXYpO1xufVxuXG4vL0NvbXB1dGVzIGxvZyBiYXNlIDIgb2YgdlxuZXhwb3J0cy5sb2cyID0gZnVuY3Rpb24odikge1xuICB2YXIgciwgc2hpZnQ7XG4gIHIgPSAgICAgKHYgPiAweEZGRkYpIDw8IDQ7IHYgPj4+PSByO1xuICBzaGlmdCA9ICh2ID4gMHhGRiAgKSA8PCAzOyB2ID4+Pj0gc2hpZnQ7IHIgfD0gc2hpZnQ7XG4gIHNoaWZ0ID0gKHYgPiAweEYgICApIDw8IDI7IHYgPj4+PSBzaGlmdDsgciB8PSBzaGlmdDtcbiAgc2hpZnQgPSAodiA+IDB4MyAgICkgPDwgMTsgdiA+Pj49IHNoaWZ0OyByIHw9IHNoaWZ0O1xuICByZXR1cm4gciB8ICh2ID4+IDEpO1xufVxuXG4vL0NvbXB1dGVzIGxvZyBiYXNlIDEwIG9mIHZcbmV4cG9ydHMubG9nMTAgPSBmdW5jdGlvbih2KSB7XG4gIHJldHVybiAgKHYgPj0gMTAwMDAwMDAwMCkgPyA5IDogKHYgPj0gMTAwMDAwMDAwKSA/IDggOiAodiA+PSAxMDAwMDAwMCkgPyA3IDpcbiAgICAgICAgICAodiA+PSAxMDAwMDAwKSA/IDYgOiAodiA+PSAxMDAwMDApID8gNSA6ICh2ID49IDEwMDAwKSA/IDQgOlxuICAgICAgICAgICh2ID49IDEwMDApID8gMyA6ICh2ID49IDEwMCkgPyAyIDogKHYgPj0gMTApID8gMSA6IDA7XG59XG5cbi8vQ291bnRzIG51bWJlciBvZiBiaXRzXG5leHBvcnRzLnBvcENvdW50ID0gZnVuY3Rpb24odikge1xuICB2ID0gdiAtICgodiA+Pj4gMSkgJiAweDU1NTU1NTU1KTtcbiAgdiA9ICh2ICYgMHgzMzMzMzMzMykgKyAoKHYgPj4+IDIpICYgMHgzMzMzMzMzMyk7XG4gIHJldHVybiAoKHYgKyAodiA+Pj4gNCkgJiAweEYwRjBGMEYpICogMHgxMDEwMTAxKSA+Pj4gMjQ7XG59XG5cbi8vQ291bnRzIG51bWJlciBvZiB0cmFpbGluZyB6ZXJvc1xuZnVuY3Rpb24gY291bnRUcmFpbGluZ1plcm9zKHYpIHtcbiAgdmFyIGMgPSAzMjtcbiAgdiAmPSAtdjtcbiAgaWYgKHYpIGMtLTtcbiAgaWYgKHYgJiAweDAwMDBGRkZGKSBjIC09IDE2O1xuICBpZiAodiAmIDB4MDBGRjAwRkYpIGMgLT0gODtcbiAgaWYgKHYgJiAweDBGMEYwRjBGKSBjIC09IDQ7XG4gIGlmICh2ICYgMHgzMzMzMzMzMykgYyAtPSAyO1xuICBpZiAodiAmIDB4NTU1NTU1NTUpIGMgLT0gMTtcbiAgcmV0dXJuIGM7XG59XG5leHBvcnRzLmNvdW50VHJhaWxpbmdaZXJvcyA9IGNvdW50VHJhaWxpbmdaZXJvcztcblxuLy9Sb3VuZHMgdG8gbmV4dCBwb3dlciBvZiAyXG5leHBvcnRzLm5leHRQb3cyID0gZnVuY3Rpb24odikge1xuICB2ICs9IHYgPT09IDA7XG4gIC0tdjtcbiAgdiB8PSB2ID4+PiAxO1xuICB2IHw9IHYgPj4+IDI7XG4gIHYgfD0gdiA+Pj4gNDtcbiAgdiB8PSB2ID4+PiA4O1xuICB2IHw9IHYgPj4+IDE2O1xuICByZXR1cm4gdiArIDE7XG59XG5cbi8vUm91bmRzIGRvd24gdG8gcHJldmlvdXMgcG93ZXIgb2YgMlxuZXhwb3J0cy5wcmV2UG93MiA9IGZ1bmN0aW9uKHYpIHtcbiAgdiB8PSB2ID4+PiAxO1xuICB2IHw9IHYgPj4+IDI7XG4gIHYgfD0gdiA+Pj4gNDtcbiAgdiB8PSB2ID4+PiA4O1xuICB2IHw9IHYgPj4+IDE2O1xuICByZXR1cm4gdiAtICh2Pj4+MSk7XG59XG5cbi8vQ29tcHV0ZXMgcGFyaXR5IG9mIHdvcmRcbmV4cG9ydHMucGFyaXR5ID0gZnVuY3Rpb24odikge1xuICB2IF49IHYgPj4+IDE2O1xuICB2IF49IHYgPj4+IDg7XG4gIHYgXj0gdiA+Pj4gNDtcbiAgdiAmPSAweGY7XG4gIHJldHVybiAoMHg2OTk2ID4+PiB2KSAmIDE7XG59XG5cbnZhciBSRVZFUlNFX1RBQkxFID0gbmV3IEFycmF5KDI1Nik7XG5cbihmdW5jdGlvbih0YWIpIHtcbiAgZm9yKHZhciBpPTA7IGk8MjU2OyArK2kpIHtcbiAgICB2YXIgdiA9IGksIHIgPSBpLCBzID0gNztcbiAgICBmb3IgKHYgPj4+PSAxOyB2OyB2ID4+Pj0gMSkge1xuICAgICAgciA8PD0gMTtcbiAgICAgIHIgfD0gdiAmIDE7XG4gICAgICAtLXM7XG4gICAgfVxuICAgIHRhYltpXSA9IChyIDw8IHMpICYgMHhmZjtcbiAgfVxufSkoUkVWRVJTRV9UQUJMRSk7XG5cbi8vUmV2ZXJzZSBiaXRzIGluIGEgMzIgYml0IHdvcmRcbmV4cG9ydHMucmV2ZXJzZSA9IGZ1bmN0aW9uKHYpIHtcbiAgcmV0dXJuICAoUkVWRVJTRV9UQUJMRVsgdiAgICAgICAgICYgMHhmZl0gPDwgMjQpIHxcbiAgICAgICAgICAoUkVWRVJTRV9UQUJMRVsodiA+Pj4gOCkgICYgMHhmZl0gPDwgMTYpIHxcbiAgICAgICAgICAoUkVWRVJTRV9UQUJMRVsodiA+Pj4gMTYpICYgMHhmZl0gPDwgOCkgIHxcbiAgICAgICAgICAgUkVWRVJTRV9UQUJMRVsodiA+Pj4gMjQpICYgMHhmZl07XG59XG5cbi8vSW50ZXJsZWF2ZSBiaXRzIG9mIDIgY29vcmRpbmF0ZXMgd2l0aCAxNiBiaXRzLiAgVXNlZnVsIGZvciBmYXN0IHF1YWR0cmVlIGNvZGVzXG5leHBvcnRzLmludGVybGVhdmUyID0gZnVuY3Rpb24oeCwgeSkge1xuICB4ICY9IDB4RkZGRjtcbiAgeCA9ICh4IHwgKHggPDwgOCkpICYgMHgwMEZGMDBGRjtcbiAgeCA9ICh4IHwgKHggPDwgNCkpICYgMHgwRjBGMEYwRjtcbiAgeCA9ICh4IHwgKHggPDwgMikpICYgMHgzMzMzMzMzMztcbiAgeCA9ICh4IHwgKHggPDwgMSkpICYgMHg1NTU1NTU1NTtcblxuICB5ICY9IDB4RkZGRjtcbiAgeSA9ICh5IHwgKHkgPDwgOCkpICYgMHgwMEZGMDBGRjtcbiAgeSA9ICh5IHwgKHkgPDwgNCkpICYgMHgwRjBGMEYwRjtcbiAgeSA9ICh5IHwgKHkgPDwgMikpICYgMHgzMzMzMzMzMztcbiAgeSA9ICh5IHwgKHkgPDwgMSkpICYgMHg1NTU1NTU1NTtcblxuICByZXR1cm4geCB8ICh5IDw8IDEpO1xufVxuXG4vL0V4dHJhY3RzIHRoZSBudGggaW50ZXJsZWF2ZWQgY29tcG9uZW50XG5leHBvcnRzLmRlaW50ZXJsZWF2ZTIgPSBmdW5jdGlvbih2LCBuKSB7XG4gIHYgPSAodiA+Pj4gbikgJiAweDU1NTU1NTU1O1xuICB2ID0gKHYgfCAodiA+Pj4gMSkpICAmIDB4MzMzMzMzMzM7XG4gIHYgPSAodiB8ICh2ID4+PiAyKSkgICYgMHgwRjBGMEYwRjtcbiAgdiA9ICh2IHwgKHYgPj4+IDQpKSAgJiAweDAwRkYwMEZGO1xuICB2ID0gKHYgfCAodiA+Pj4gMTYpKSAmIDB4MDAwRkZGRjtcbiAgcmV0dXJuICh2IDw8IDE2KSA+PiAxNjtcbn1cblxuXG4vL0ludGVybGVhdmUgYml0cyBvZiAzIGNvb3JkaW5hdGVzLCBlYWNoIHdpdGggMTAgYml0cy4gIFVzZWZ1bCBmb3IgZmFzdCBvY3RyZWUgY29kZXNcbmV4cG9ydHMuaW50ZXJsZWF2ZTMgPSBmdW5jdGlvbih4LCB5LCB6KSB7XG4gIHggJj0gMHgzRkY7XG4gIHggID0gKHggfCAoeDw8MTYpKSAmIDQyNzgxOTAzMzU7XG4gIHggID0gKHggfCAoeDw8OCkpICAmIDI1MTcxOTY5NTtcbiAgeCAgPSAoeCB8ICh4PDw0KSkgICYgMzI3MjM1NjAzNTtcbiAgeCAgPSAoeCB8ICh4PDwyKSkgICYgMTIyNzEzMzUxMztcblxuICB5ICY9IDB4M0ZGO1xuICB5ICA9ICh5IHwgKHk8PDE2KSkgJiA0Mjc4MTkwMzM1O1xuICB5ICA9ICh5IHwgKHk8PDgpKSAgJiAyNTE3MTk2OTU7XG4gIHkgID0gKHkgfCAoeTw8NCkpICAmIDMyNzIzNTYwMzU7XG4gIHkgID0gKHkgfCAoeTw8MikpICAmIDEyMjcxMzM1MTM7XG4gIHggfD0gKHkgPDwgMSk7XG4gIFxuICB6ICY9IDB4M0ZGO1xuICB6ICA9ICh6IHwgKHo8PDE2KSkgJiA0Mjc4MTkwMzM1O1xuICB6ICA9ICh6IHwgKHo8PDgpKSAgJiAyNTE3MTk2OTU7XG4gIHogID0gKHogfCAoejw8NCkpICAmIDMyNzIzNTYwMzU7XG4gIHogID0gKHogfCAoejw8MikpICAmIDEyMjcxMzM1MTM7XG4gIFxuICByZXR1cm4geCB8ICh6IDw8IDIpO1xufVxuXG4vL0V4dHJhY3RzIG50aCBpbnRlcmxlYXZlZCBjb21wb25lbnQgb2YgYSAzLXR1cGxlXG5leHBvcnRzLmRlaW50ZXJsZWF2ZTMgPSBmdW5jdGlvbih2LCBuKSB7XG4gIHYgPSAodiA+Pj4gbikgICAgICAgJiAxMjI3MTMzNTEzO1xuICB2ID0gKHYgfCAodj4+PjIpKSAgICYgMzI3MjM1NjAzNTtcbiAgdiA9ICh2IHwgKHY+Pj40KSkgICAmIDI1MTcxOTY5NTtcbiAgdiA9ICh2IHwgKHY+Pj44KSkgICAmIDQyNzgxOTAzMzU7XG4gIHYgPSAodiB8ICh2Pj4+MTYpKSAgJiAweDNGRjtcbiAgcmV0dXJuICh2PDwyMik+PjIyO1xufVxuXG4vL0NvbXB1dGVzIG5leHQgY29tYmluYXRpb24gaW4gY29sZXhpY29ncmFwaGljIG9yZGVyICh0aGlzIGlzIG1pc3Rha2VubHkgY2FsbGVkIG5leHRQZXJtdXRhdGlvbiBvbiB0aGUgYml0IHR3aWRkbGluZyBoYWNrcyBwYWdlKVxuZXhwb3J0cy5uZXh0Q29tYmluYXRpb24gPSBmdW5jdGlvbih2KSB7XG4gIHZhciB0ID0gdiB8ICh2IC0gMSk7XG4gIHJldHVybiAodCArIDEpIHwgKCgofnQgJiAtfnQpIC0gMSkgPj4+IChjb3VudFRyYWlsaW5nWmVyb3ModikgKyAxKSk7XG59XG5cbiIsIlwidXNlIHN0cmljdFwiOyBcInVzZSByZXN0cmljdFwiO1xuXG5tb2R1bGUuZXhwb3J0cyA9IFVuaW9uRmluZDtcblxuZnVuY3Rpb24gVW5pb25GaW5kKGNvdW50KSB7XG4gIHRoaXMucm9vdHMgPSBuZXcgQXJyYXkoY291bnQpO1xuICB0aGlzLnJhbmtzID0gbmV3IEFycmF5KGNvdW50KTtcbiAgXG4gIGZvcih2YXIgaT0wOyBpPGNvdW50OyArK2kpIHtcbiAgICB0aGlzLnJvb3RzW2ldID0gaTtcbiAgICB0aGlzLnJhbmtzW2ldID0gMDtcbiAgfVxufVxuXG52YXIgcHJvdG8gPSBVbmlvbkZpbmQucHJvdG90eXBlXG5cbk9iamVjdC5kZWZpbmVQcm9wZXJ0eShwcm90bywgXCJsZW5ndGhcIiwge1xuICBcImdldFwiOiBmdW5jdGlvbigpIHtcbiAgICByZXR1cm4gdGhpcy5yb290cy5sZW5ndGhcbiAgfVxufSlcblxucHJvdG8ubWFrZVNldCA9IGZ1bmN0aW9uKCkge1xuICB2YXIgbiA9IHRoaXMucm9vdHMubGVuZ3RoO1xuICB0aGlzLnJvb3RzLnB1c2gobik7XG4gIHRoaXMucmFua3MucHVzaCgwKTtcbiAgcmV0dXJuIG47XG59XG5cbnByb3RvLmZpbmQgPSBmdW5jdGlvbih4KSB7XG4gIHZhciByb290cyA9IHRoaXMucm9vdHM7XG4gIHdoaWxlKHJvb3RzW3hdICE9PSB4KSB7XG4gICAgdmFyIHkgPSByb290c1t4XTtcbiAgICByb290c1t4XSA9IHJvb3RzW3ldO1xuICAgIHggPSB5O1xuICB9XG4gIHJldHVybiB4O1xufVxuXG5wcm90by5saW5rID0gZnVuY3Rpb24oeCwgeSkge1xuICB2YXIgeHIgPSB0aGlzLmZpbmQoeClcbiAgICAsIHlyID0gdGhpcy5maW5kKHkpO1xuICBpZih4ciA9PT0geXIpIHtcbiAgICByZXR1cm47XG4gIH1cbiAgdmFyIHJhbmtzID0gdGhpcy5yYW5rc1xuICAgICwgcm9vdHMgPSB0aGlzLnJvb3RzXG4gICAgLCB4ZCAgICA9IHJhbmtzW3hyXVxuICAgICwgeWQgICAgPSByYW5rc1t5cl07XG4gIGlmKHhkIDwgeWQpIHtcbiAgICByb290c1t4cl0gPSB5cjtcbiAgfSBlbHNlIGlmKHlkIDwgeGQpIHtcbiAgICByb290c1t5cl0gPSB4cjtcbiAgfSBlbHNlIHtcbiAgICByb290c1t5cl0gPSB4cjtcbiAgICArK3JhbmtzW3hyXTtcbiAgfVxufSIsIlwidXNlIHN0cmljdFwiOyBcInVzZSByZXN0cmljdFwiO1xuXG52YXIgYml0cyAgICAgID0gcmVxdWlyZShcImJpdC10d2lkZGxlXCIpXG4gICwgVW5pb25GaW5kID0gcmVxdWlyZShcInVuaW9uLWZpbmRcIilcblxuLy9SZXR1cm5zIHRoZSBkaW1lbnNpb24gb2YgYSBjZWxsIGNvbXBsZXhcbmZ1bmN0aW9uIGRpbWVuc2lvbihjZWxscykge1xuICB2YXIgZCA9IDBcbiAgICAsIG1heCA9IE1hdGgubWF4XG4gIGZvcih2YXIgaT0wLCBpbD1jZWxscy5sZW5ndGg7IGk8aWw7ICsraSkge1xuICAgIGQgPSBtYXgoZCwgY2VsbHNbaV0ubGVuZ3RoKVxuICB9XG4gIHJldHVybiBkLTFcbn1cbmV4cG9ydHMuZGltZW5zaW9uID0gZGltZW5zaW9uXG5cbi8vQ291bnRzIHRoZSBudW1iZXIgb2YgdmVydGljZXMgaW4gZmFjZXNcbmZ1bmN0aW9uIGNvdW50VmVydGljZXMoY2VsbHMpIHtcbiAgdmFyIHZjID0gLTFcbiAgICAsIG1heCA9IE1hdGgubWF4XG4gIGZvcih2YXIgaT0wLCBpbD1jZWxscy5sZW5ndGg7IGk8aWw7ICsraSkge1xuICAgIHZhciBjID0gY2VsbHNbaV1cbiAgICBmb3IodmFyIGo9MCwgamw9Yy5sZW5ndGg7IGo8amw7ICsraikge1xuICAgICAgdmMgPSBtYXgodmMsIGNbal0pXG4gICAgfVxuICB9XG4gIHJldHVybiB2YysxXG59XG5leHBvcnRzLmNvdW50VmVydGljZXMgPSBjb3VudFZlcnRpY2VzXG5cbi8vUmV0dXJucyBhIGRlZXAgY29weSBvZiBjZWxsc1xuZnVuY3Rpb24gY2xvbmVDZWxscyhjZWxscykge1xuICB2YXIgbmNlbGxzID0gbmV3IEFycmF5KGNlbGxzLmxlbmd0aClcbiAgZm9yKHZhciBpPTAsIGlsPWNlbGxzLmxlbmd0aDsgaTxpbDsgKytpKSB7XG4gICAgbmNlbGxzW2ldID0gY2VsbHNbaV0uc2xpY2UoMClcbiAgfVxuICByZXR1cm4gbmNlbGxzXG59XG5leHBvcnRzLmNsb25lQ2VsbHMgPSBjbG9uZUNlbGxzXG5cbi8vUmFua3MgYSBwYWlyIG9mIGNlbGxzIHVwIHRvIHBlcm11dGF0aW9uXG5mdW5jdGlvbiBjb21wYXJlQ2VsbHMoYSwgYikge1xuICB2YXIgbiA9IGEubGVuZ3RoXG4gICAgLCB0ID0gYS5sZW5ndGggLSBiLmxlbmd0aFxuICAgICwgbWluID0gTWF0aC5taW5cbiAgaWYodCkge1xuICAgIHJldHVybiB0XG4gIH1cbiAgc3dpdGNoKG4pIHtcbiAgICBjYXNlIDA6XG4gICAgICByZXR1cm4gMDtcbiAgICBjYXNlIDE6XG4gICAgICByZXR1cm4gYVswXSAtIGJbMF07XG4gICAgY2FzZSAyOlxuICAgICAgdmFyIGQgPSBhWzBdK2FbMV0tYlswXS1iWzFdXG4gICAgICBpZihkKSB7XG4gICAgICAgIHJldHVybiBkXG4gICAgICB9XG4gICAgICByZXR1cm4gbWluKGFbMF0sYVsxXSkgLSBtaW4oYlswXSxiWzFdKVxuICAgIGNhc2UgMzpcbiAgICAgIHZhciBsMSA9IGFbMF0rYVsxXVxuICAgICAgICAsIG0xID0gYlswXStiWzFdXG4gICAgICBkID0gbDErYVsyXSAtIChtMStiWzJdKVxuICAgICAgaWYoZCkge1xuICAgICAgICByZXR1cm4gZFxuICAgICAgfVxuICAgICAgdmFyIGwwID0gbWluKGFbMF0sIGFbMV0pXG4gICAgICAgICwgbTAgPSBtaW4oYlswXSwgYlsxXSlcbiAgICAgICAgLCBkICA9IG1pbihsMCwgYVsyXSkgLSBtaW4obTAsIGJbMl0pXG4gICAgICBpZihkKSB7XG4gICAgICAgIHJldHVybiBkXG4gICAgICB9XG4gICAgICByZXR1cm4gbWluKGwwK2FbMl0sIGwxKSAtIG1pbihtMCtiWzJdLCBtMSlcbiAgICBcbiAgICAvL1RPRE86IE1heWJlIG9wdGltaXplIG49NCBhcyB3ZWxsP1xuICAgIFxuICAgIGRlZmF1bHQ6XG4gICAgICB2YXIgYXMgPSBhLnNsaWNlKDApXG4gICAgICBhcy5zb3J0KClcbiAgICAgIHZhciBicyA9IGIuc2xpY2UoMClcbiAgICAgIGJzLnNvcnQoKVxuICAgICAgZm9yKHZhciBpPTA7IGk8bjsgKytpKSB7XG4gICAgICAgIHQgPSBhc1tpXSAtIGJzW2ldXG4gICAgICAgIGlmKHQpIHtcbiAgICAgICAgICByZXR1cm4gdFxuICAgICAgICB9XG4gICAgICB9XG4gICAgICByZXR1cm4gMFxuICB9XG59XG5leHBvcnRzLmNvbXBhcmVDZWxscyA9IGNvbXBhcmVDZWxsc1xuXG5mdW5jdGlvbiBjb21wYXJlWmlwcGVkKGEsIGIpIHtcbiAgcmV0dXJuIGNvbXBhcmVDZWxscyhhWzBdLCBiWzBdKVxufVxuXG4vL1B1dHMgYSBjZWxsIGNvbXBsZXggaW50byBub3JtYWwgb3JkZXIgZm9yIHRoZSBwdXJwb3NlcyBvZiBmaW5kQ2VsbCBxdWVyaWVzXG5mdW5jdGlvbiBub3JtYWxpemUoY2VsbHMsIGF0dHIpIHtcbiAgaWYoYXR0cikge1xuICAgIHZhciBsZW4gPSBjZWxscy5sZW5ndGhcbiAgICB2YXIgemlwcGVkID0gbmV3IEFycmF5KGxlbilcbiAgICBmb3IodmFyIGk9MDsgaTxsZW47ICsraSkge1xuICAgICAgemlwcGVkW2ldID0gW2NlbGxzW2ldLCBhdHRyW2ldXVxuICAgIH1cbiAgICB6aXBwZWQuc29ydChjb21wYXJlWmlwcGVkKVxuICAgIGZvcih2YXIgaT0wOyBpPGxlbjsgKytpKSB7XG4gICAgICBjZWxsc1tpXSA9IHppcHBlZFtpXVswXVxuICAgICAgYXR0cltpXSA9IHppcHBlZFtpXVsxXVxuICAgIH1cbiAgICByZXR1cm4gY2VsbHNcbiAgfSBlbHNlIHtcbiAgICBjZWxscy5zb3J0KGNvbXBhcmVDZWxscylcbiAgICByZXR1cm4gY2VsbHNcbiAgfVxufVxuZXhwb3J0cy5ub3JtYWxpemUgPSBub3JtYWxpemVcblxuLy9SZW1vdmVzIGFsbCBkdXBsaWNhdGUgY2VsbHMgaW4gdGhlIGNvbXBsZXhcbmZ1bmN0aW9uIHVuaXF1ZShjZWxscykge1xuICBpZihjZWxscy5sZW5ndGggPT09IDApIHtcbiAgICByZXR1cm4gW11cbiAgfVxuICB2YXIgcHRyID0gMVxuICAgICwgbGVuID0gY2VsbHMubGVuZ3RoXG4gIGZvcih2YXIgaT0xOyBpPGxlbjsgKytpKSB7XG4gICAgdmFyIGEgPSBjZWxsc1tpXVxuICAgIGlmKGNvbXBhcmVDZWxscyhhLCBjZWxsc1tpLTFdKSkge1xuICAgICAgaWYoaSA9PT0gcHRyKSB7XG4gICAgICAgIHB0cisrXG4gICAgICAgIGNvbnRpbnVlXG4gICAgICB9XG4gICAgICBjZWxsc1twdHIrK10gPSBhXG4gICAgfVxuICB9XG4gIGNlbGxzLmxlbmd0aCA9IHB0clxuICByZXR1cm4gY2VsbHNcbn1cbmV4cG9ydHMudW5pcXVlID0gdW5pcXVlO1xuXG4vL0ZpbmRzIGEgY2VsbCBpbiBhIG5vcm1hbGl6ZWQgY2VsbCBjb21wbGV4XG5mdW5jdGlvbiBmaW5kQ2VsbChjZWxscywgYykge1xuICB2YXIgbG8gPSAwXG4gICAgLCBoaSA9IGNlbGxzLmxlbmd0aC0xXG4gICAgLCByICA9IC0xXG4gIHdoaWxlIChsbyA8PSBoaSkge1xuICAgIHZhciBtaWQgPSAobG8gKyBoaSkgPj4gMVxuICAgICAgLCBzICAgPSBjb21wYXJlQ2VsbHMoY2VsbHNbbWlkXSwgYylcbiAgICBpZihzIDw9IDApIHtcbiAgICAgIGlmKHMgPT09IDApIHtcbiAgICAgICAgciA9IG1pZFxuICAgICAgfVxuICAgICAgbG8gPSBtaWQgKyAxXG4gICAgfSBlbHNlIGlmKHMgPiAwKSB7XG4gICAgICBoaSA9IG1pZCAtIDFcbiAgICB9XG4gIH1cbiAgcmV0dXJuIHJcbn1cbmV4cG9ydHMuZmluZENlbGwgPSBmaW5kQ2VsbDtcblxuLy9CdWlsZHMgYW4gaW5kZXggZm9yIGFuIG4tY2VsbC4gIFRoaXMgaXMgbW9yZSBnZW5lcmFsIHRoYW4gZHVhbCwgYnV0IGxlc3MgZWZmaWNpZW50XG5mdW5jdGlvbiBpbmNpZGVuY2UoZnJvbV9jZWxscywgdG9fY2VsbHMpIHtcbiAgdmFyIGluZGV4ID0gbmV3IEFycmF5KGZyb21fY2VsbHMubGVuZ3RoKVxuICBmb3IodmFyIGk9MCwgaWw9aW5kZXgubGVuZ3RoOyBpPGlsOyArK2kpIHtcbiAgICBpbmRleFtpXSA9IFtdXG4gIH1cbiAgdmFyIGIgPSBbXVxuICBmb3IodmFyIGk9MCwgbj10b19jZWxscy5sZW5ndGg7IGk8bjsgKytpKSB7XG4gICAgdmFyIGMgPSB0b19jZWxsc1tpXVxuICAgIHZhciBjbCA9IGMubGVuZ3RoXG4gICAgZm9yKHZhciBrPTEsIGtuPSgxPDxjbCk7IGs8a247ICsraykge1xuICAgICAgYi5sZW5ndGggPSBiaXRzLnBvcENvdW50KGspXG4gICAgICB2YXIgbCA9IDBcbiAgICAgIGZvcih2YXIgaj0wOyBqPGNsOyArK2opIHtcbiAgICAgICAgaWYoayAmICgxPDxqKSkge1xuICAgICAgICAgIGJbbCsrXSA9IGNbal1cbiAgICAgICAgfVxuICAgICAgfVxuICAgICAgdmFyIGlkeD1maW5kQ2VsbChmcm9tX2NlbGxzLCBiKVxuICAgICAgaWYoaWR4IDwgMCkge1xuICAgICAgICBjb250aW51ZVxuICAgICAgfVxuICAgICAgd2hpbGUodHJ1ZSkge1xuICAgICAgICBpbmRleFtpZHgrK10ucHVzaChpKVxuICAgICAgICBpZihpZHggPj0gZnJvbV9jZWxscy5sZW5ndGggfHwgY29tcGFyZUNlbGxzKGZyb21fY2VsbHNbaWR4XSwgYikgIT09IDApIHtcbiAgICAgICAgICBicmVha1xuICAgICAgICB9XG4gICAgICB9XG4gICAgfVxuICB9XG4gIHJldHVybiBpbmRleFxufVxuZXhwb3J0cy5pbmNpZGVuY2UgPSBpbmNpZGVuY2VcblxuLy9Db21wdXRlcyB0aGUgZHVhbCBvZiB0aGUgbWVzaC4gIFRoaXMgaXMgYmFzaWNhbGx5IGFuIG9wdGltaXplZCB2ZXJzaW9uIG9mIGJ1aWxkSW5kZXggZm9yIHRoZSBzaXR1YXRpb24gd2hlcmUgZnJvbV9jZWxscyBpcyBqdXN0IHRoZSBsaXN0IG9mIHZlcnRpY2VzXG5mdW5jdGlvbiBkdWFsKGNlbGxzLCB2ZXJ0ZXhfY291bnQpIHtcbiAgaWYoIXZlcnRleF9jb3VudCkge1xuICAgIHJldHVybiBpbmNpZGVuY2UodW5pcXVlKHNrZWxldG9uKGNlbGxzLCAwKSksIGNlbGxzLCAwKVxuICB9XG4gIHZhciByZXMgPSBuZXcgQXJyYXkodmVydGV4X2NvdW50KVxuICBmb3IodmFyIGk9MDsgaTx2ZXJ0ZXhfY291bnQ7ICsraSkge1xuICAgIHJlc1tpXSA9IFtdXG4gIH1cbiAgZm9yKHZhciBpPTAsIGxlbj1jZWxscy5sZW5ndGg7IGk8bGVuOyArK2kpIHtcbiAgICB2YXIgYyA9IGNlbGxzW2ldXG4gICAgZm9yKHZhciBqPTAsIGNsPWMubGVuZ3RoOyBqPGNsOyArK2opIHtcbiAgICAgIHJlc1tjW2pdXS5wdXNoKGkpXG4gICAgfVxuICB9XG4gIHJldHVybiByZXNcbn1cbmV4cG9ydHMuZHVhbCA9IGR1YWxcblxuLy9FbnVtZXJhdGVzIGFsbCBjZWxscyBpbiB0aGUgY29tcGxleFxuZnVuY3Rpb24gZXhwbG9kZShjZWxscykge1xuICB2YXIgcmVzdWx0ID0gW11cbiAgZm9yKHZhciBpPTAsIGlsPWNlbGxzLmxlbmd0aDsgaTxpbDsgKytpKSB7XG4gICAgdmFyIGMgPSBjZWxsc1tpXVxuICAgICAgLCBjbCA9IGMubGVuZ3RofDBcbiAgICBmb3IodmFyIGo9MSwgamw9KDE8PGNsKTsgajxqbDsgKytqKSB7XG4gICAgICB2YXIgYiA9IFtdXG4gICAgICBmb3IodmFyIGs9MDsgazxjbDsgKytrKSB7XG4gICAgICAgIGlmKChqID4+PiBrKSAmIDEpIHtcbiAgICAgICAgICBiLnB1c2goY1trXSlcbiAgICAgICAgfVxuICAgICAgfVxuICAgICAgcmVzdWx0LnB1c2goYilcbiAgICB9XG4gIH1cbiAgcmV0dXJuIG5vcm1hbGl6ZShyZXN1bHQpXG59XG5leHBvcnRzLmV4cGxvZGUgPSBleHBsb2RlXG5cbi8vRW51bWVyYXRlcyBhbGwgb2YgdGhlIG4tY2VsbHMgb2YgYSBjZWxsIGNvbXBsZXhcbmZ1bmN0aW9uIHNrZWxldG9uKGNlbGxzLCBuKSB7XG4gIGlmKG4gPCAwKSB7XG4gICAgcmV0dXJuIFtdXG4gIH1cbiAgdmFyIHJlc3VsdCA9IFtdXG4gICAgLCBrMCAgICAgPSAoMTw8KG4rMSkpLTFcbiAgZm9yKHZhciBpPTA7IGk8Y2VsbHMubGVuZ3RoOyArK2kpIHtcbiAgICB2YXIgYyA9IGNlbGxzW2ldXG4gICAgZm9yKHZhciBrPWswOyBrPCgxPDxjLmxlbmd0aCk7IGs9Yml0cy5uZXh0Q29tYmluYXRpb24oaykpIHtcbiAgICAgIHZhciBiID0gbmV3IEFycmF5KG4rMSlcbiAgICAgICAgLCBsID0gMFxuICAgICAgZm9yKHZhciBqPTA7IGo8Yy5sZW5ndGg7ICsraikge1xuICAgICAgICBpZihrICYgKDE8PGopKSB7XG4gICAgICAgICAgYltsKytdID0gY1tqXVxuICAgICAgICB9XG4gICAgICB9XG4gICAgICByZXN1bHQucHVzaChiKVxuICAgIH1cbiAgfVxuICByZXR1cm4gbm9ybWFsaXplKHJlc3VsdClcbn1cbmV4cG9ydHMuc2tlbGV0b24gPSBza2VsZXRvbjtcblxuLy9Db21wdXRlcyB0aGUgYm91bmRhcnkgb2YgYWxsIGNlbGxzLCBkb2VzIG5vdCByZW1vdmUgZHVwbGljYXRlc1xuZnVuY3Rpb24gYm91bmRhcnkoY2VsbHMpIHtcbiAgdmFyIHJlcyA9IFtdXG4gIGZvcih2YXIgaT0wLGlsPWNlbGxzLmxlbmd0aDsgaTxpbDsgKytpKSB7XG4gICAgdmFyIGMgPSBjZWxsc1tpXVxuICAgIGZvcih2YXIgaj0wLGNsPWMubGVuZ3RoOyBqPGNsOyArK2opIHtcbiAgICAgIHZhciBiID0gbmV3IEFycmF5KGMubGVuZ3RoLTEpXG4gICAgICBmb3IodmFyIGs9MCwgbD0wOyBrPGNsOyArK2spIHtcbiAgICAgICAgaWYoayAhPT0gaikge1xuICAgICAgICAgIGJbbCsrXSA9IGNba11cbiAgICAgICAgfVxuICAgICAgfVxuICAgICAgcmVzLnB1c2goYilcbiAgICB9XG4gIH1cbiAgcmV0dXJuIG5vcm1hbGl6ZShyZXMpXG59XG5leHBvcnRzLmJvdW5kYXJ5ID0gYm91bmRhcnk7XG5cbi8vQ29tcHV0ZXMgY29ubmVjdGVkIGNvbXBvbmVudHMgZm9yIGEgZGVuc2UgY2VsbCBjb21wbGV4XG5mdW5jdGlvbiBjb25uZWN0ZWRDb21wb25lbnRzX2RlbnNlKGNlbGxzLCB2ZXJ0ZXhfY291bnQpIHtcbiAgdmFyIGxhYmVscyA9IG5ldyBVbmlvbkZpbmQodmVydGV4X2NvdW50KVxuICBmb3IodmFyIGk9MDsgaTxjZWxscy5sZW5ndGg7ICsraSkge1xuICAgIHZhciBjID0gY2VsbHNbaV1cbiAgICBmb3IodmFyIGo9MDsgajxjLmxlbmd0aDsgKytqKSB7XG4gICAgICBmb3IodmFyIGs9aisxOyBrPGMubGVuZ3RoOyArK2spIHtcbiAgICAgICAgbGFiZWxzLmxpbmsoY1tqXSwgY1trXSlcbiAgICAgIH1cbiAgICB9XG4gIH1cbiAgdmFyIGNvbXBvbmVudHMgPSBbXVxuICAgICwgY29tcG9uZW50X2xhYmVscyA9IGxhYmVscy5yYW5rc1xuICBmb3IodmFyIGk9MDsgaTxjb21wb25lbnRfbGFiZWxzLmxlbmd0aDsgKytpKSB7XG4gICAgY29tcG9uZW50X2xhYmVsc1tpXSA9IC0xXG4gIH1cbiAgZm9yKHZhciBpPTA7IGk8Y2VsbHMubGVuZ3RoOyArK2kpIHtcbiAgICB2YXIgbCA9IGxhYmVscy5maW5kKGNlbGxzW2ldWzBdKVxuICAgIGlmKGNvbXBvbmVudF9sYWJlbHNbbF0gPCAwKSB7XG4gICAgICBjb21wb25lbnRfbGFiZWxzW2xdID0gY29tcG9uZW50cy5sZW5ndGhcbiAgICAgIGNvbXBvbmVudHMucHVzaChbY2VsbHNbaV0uc2xpY2UoMCldKVxuICAgIH0gZWxzZSB7XG4gICAgICBjb21wb25lbnRzW2NvbXBvbmVudF9sYWJlbHNbbF1dLnB1c2goY2VsbHNbaV0uc2xpY2UoMCkpXG4gICAgfVxuICB9XG4gIHJldHVybiBjb21wb25lbnRzXG59XG5cbi8vQ29tcHV0ZXMgY29ubmVjdGVkIGNvbXBvbmVudHMgZm9yIGEgc3BhcnNlIGdyYXBoXG5mdW5jdGlvbiBjb25uZWN0ZWRDb21wb25lbnRzX3NwYXJzZShjZWxscykge1xuICB2YXIgdmVydGljZXMgID0gdW5pcXVlKG5vcm1hbGl6ZShza2VsZXRvbihjZWxscywgMCkpKVxuICAgICwgbGFiZWxzICAgID0gbmV3IFVuaW9uRmluZCh2ZXJ0aWNlcy5sZW5ndGgpXG4gIGZvcih2YXIgaT0wOyBpPGNlbGxzLmxlbmd0aDsgKytpKSB7XG4gICAgdmFyIGMgPSBjZWxsc1tpXVxuICAgIGZvcih2YXIgaj0wOyBqPGMubGVuZ3RoOyArK2opIHtcbiAgICAgIHZhciB2aiA9IGZpbmRDZWxsKHZlcnRpY2VzLCBbY1tqXV0pXG4gICAgICBmb3IodmFyIGs9aisxOyBrPGMubGVuZ3RoOyArK2spIHtcbiAgICAgICAgbGFiZWxzLmxpbmsodmosIGZpbmRDZWxsKHZlcnRpY2VzLCBbY1trXV0pKVxuICAgICAgfVxuICAgIH1cbiAgfVxuICB2YXIgY29tcG9uZW50cyAgICAgICAgPSBbXVxuICAgICwgY29tcG9uZW50X2xhYmVscyAgPSBsYWJlbHMucmFua3NcbiAgZm9yKHZhciBpPTA7IGk8Y29tcG9uZW50X2xhYmVscy5sZW5ndGg7ICsraSkge1xuICAgIGNvbXBvbmVudF9sYWJlbHNbaV0gPSAtMVxuICB9XG4gIGZvcih2YXIgaT0wOyBpPGNlbGxzLmxlbmd0aDsgKytpKSB7XG4gICAgdmFyIGwgPSBsYWJlbHMuZmluZChmaW5kQ2VsbCh2ZXJ0aWNlcywgW2NlbGxzW2ldWzBdXSkpO1xuICAgIGlmKGNvbXBvbmVudF9sYWJlbHNbbF0gPCAwKSB7XG4gICAgICBjb21wb25lbnRfbGFiZWxzW2xdID0gY29tcG9uZW50cy5sZW5ndGhcbiAgICAgIGNvbXBvbmVudHMucHVzaChbY2VsbHNbaV0uc2xpY2UoMCldKVxuICAgIH0gZWxzZSB7XG4gICAgICBjb21wb25lbnRzW2NvbXBvbmVudF9sYWJlbHNbbF1dLnB1c2goY2VsbHNbaV0uc2xpY2UoMCkpXG4gICAgfVxuICB9XG4gIHJldHVybiBjb21wb25lbnRzXG59XG5cbi8vQ29tcHV0ZXMgY29ubmVjdGVkIGNvbXBvbmVudHMgZm9yIGEgY2VsbCBjb21wbGV4XG5mdW5jdGlvbiBjb25uZWN0ZWRDb21wb25lbnRzKGNlbGxzLCB2ZXJ0ZXhfY291bnQpIHtcbiAgaWYodmVydGV4X2NvdW50KSB7XG4gICAgcmV0dXJuIGNvbm5lY3RlZENvbXBvbmVudHNfZGVuc2UoY2VsbHMsIHZlcnRleF9jb3VudClcbiAgfVxuICByZXR1cm4gY29ubmVjdGVkQ29tcG9uZW50c19zcGFyc2UoY2VsbHMpXG59XG5leHBvcnRzLmNvbm5lY3RlZENvbXBvbmVudHMgPSBjb25uZWN0ZWRDb21wb25lbnRzXG4iLCJcInVzZSBzdHJpY3RcIlxuXG5mdW5jdGlvbiB1bmlxdWVfcHJlZChsaXN0LCBjb21wYXJlKSB7XG4gIHZhciBwdHIgPSAxXG4gICAgLCBsZW4gPSBsaXN0Lmxlbmd0aFxuICAgICwgYT1saXN0WzBdLCBiPWxpc3RbMF1cbiAgZm9yKHZhciBpPTE7IGk8bGVuOyArK2kpIHtcbiAgICBiID0gYVxuICAgIGEgPSBsaXN0W2ldXG4gICAgaWYoY29tcGFyZShhLCBiKSkge1xuICAgICAgaWYoaSA9PT0gcHRyKSB7XG4gICAgICAgIHB0cisrXG4gICAgICAgIGNvbnRpbnVlXG4gICAgICB9XG4gICAgICBsaXN0W3B0cisrXSA9IGFcbiAgICB9XG4gIH1cbiAgbGlzdC5sZW5ndGggPSBwdHJcbiAgcmV0dXJuIGxpc3Rcbn1cblxuZnVuY3Rpb24gdW5pcXVlX2VxKGxpc3QpIHtcbiAgdmFyIHB0ciA9IDFcbiAgICAsIGxlbiA9IGxpc3QubGVuZ3RoXG4gICAgLCBhPWxpc3RbMF0sIGIgPSBsaXN0WzBdXG4gIGZvcih2YXIgaT0xOyBpPGxlbjsgKytpLCBiPWEpIHtcbiAgICBiID0gYVxuICAgIGEgPSBsaXN0W2ldXG4gICAgaWYoYSAhPT0gYikge1xuICAgICAgaWYoaSA9PT0gcHRyKSB7XG4gICAgICAgIHB0cisrXG4gICAgICAgIGNvbnRpbnVlXG4gICAgICB9XG4gICAgICBsaXN0W3B0cisrXSA9IGFcbiAgICB9XG4gIH1cbiAgbGlzdC5sZW5ndGggPSBwdHJcbiAgcmV0dXJuIGxpc3Rcbn1cblxuZnVuY3Rpb24gdW5pcXVlKGxpc3QsIGNvbXBhcmUsIHNvcnRlZCkge1xuICBpZihsaXN0Lmxlbmd0aCA9PT0gMCkge1xuICAgIHJldHVybiBsaXN0XG4gIH1cbiAgaWYoY29tcGFyZSkge1xuICAgIGlmKCFzb3J0ZWQpIHtcbiAgICAgIGxpc3Quc29ydChjb21wYXJlKVxuICAgIH1cbiAgICByZXR1cm4gdW5pcXVlX3ByZWQobGlzdCwgY29tcGFyZSlcbiAgfVxuICBpZighc29ydGVkKSB7XG4gICAgbGlzdC5zb3J0KClcbiAgfVxuICByZXR1cm4gdW5pcXVlX2VxKGxpc3QpXG59XG5cbm1vZHVsZS5leHBvcnRzID0gdW5pcXVlXG4iLCJcInVzZSBzdHJpY3RcIlxuXG52YXIgY2ggPSByZXF1aXJlKFwiaW5jcmVtZW50YWwtY29udmV4LWh1bGxcIilcbnZhciB1bmlxID0gcmVxdWlyZShcInVuaXFcIilcblxubW9kdWxlLmV4cG9ydHMgPSB0cmlhbmd1bGF0ZVxuXG5mdW5jdGlvbiBMaWZ0ZWRQb2ludChwLCBpKSB7XG4gIHRoaXMucG9pbnQgPSBwXG4gIHRoaXMuaW5kZXggPSBpXG59XG5cbmZ1bmN0aW9uIGNvbXBhcmVMaWZ0ZWQoYSwgYikge1xuICB2YXIgYXAgPSBhLnBvaW50XG4gIHZhciBicCA9IGIucG9pbnRcbiAgdmFyIGQgPSBhcC5sZW5ndGhcbiAgZm9yKHZhciBpPTA7IGk8ZDsgKytpKSB7XG4gICAgdmFyIHMgPSBicFtpXSAtIGFwW2ldXG4gICAgaWYocykge1xuICAgICAgcmV0dXJuIHNcbiAgICB9XG4gIH1cbiAgcmV0dXJuIDBcbn1cblxuZnVuY3Rpb24gdHJpYW5ndWxhdGUxRChuLCBwb2ludHMsIGluY2x1ZGVQb2ludEF0SW5maW5pdHkpIHtcbiAgaWYobiA9PT0gMSkge1xuICAgIGlmKGluY2x1ZGVQb2ludEF0SW5maW5pdHkpIHtcbiAgICAgIHJldHVybiBbIFstMSwgMF0gXVxuICAgIH0gZWxzZSB7XG4gICAgICByZXR1cm4gW11cbiAgICB9XG4gIH1cbiAgdmFyIGxpZnRlZCA9IHBvaW50cy5tYXAoZnVuY3Rpb24ocCwgaSkge1xuICAgIHJldHVybiBbIHBbMF0sIGkgXVxuICB9KVxuICBsaWZ0ZWQuc29ydChmdW5jdGlvbihhLGIpIHtcbiAgICByZXR1cm4gYVswXSAtIGJbMF1cbiAgfSlcbiAgdmFyIGNlbGxzID0gbmV3IEFycmF5KG4gLSAxKVxuICBmb3IodmFyIGk9MTsgaTxuOyArK2kpIHtcbiAgICB2YXIgYSA9IGxpZnRlZFtpLTFdXG4gICAgdmFyIGIgPSBsaWZ0ZWRbaV1cbiAgICBjZWxsc1tpLTFdID0gWyBhWzFdLCBiWzFdIF1cbiAgfVxuICBpZihpbmNsdWRlUG9pbnRBdEluZmluaXR5KSB7XG4gICAgY2VsbHMucHVzaChcbiAgICAgIFsgLTEsIGNlbGxzWzBdWzFdLCBdLFxuICAgICAgWyBjZWxsc1tuLTFdWzFdLCAtMSBdKVxuICB9XG4gIHJldHVybiBjZWxsc1xufVxuXG5mdW5jdGlvbiB0cmlhbmd1bGF0ZShwb2ludHMsIGluY2x1ZGVQb2ludEF0SW5maW5pdHkpIHtcbiAgdmFyIG4gPSBwb2ludHMubGVuZ3RoXG4gIGlmKG4gPT09IDApIHtcbiAgICByZXR1cm4gW11cbiAgfVxuICBcbiAgdmFyIGQgPSBwb2ludHNbMF0ubGVuZ3RoXG4gIGlmKGQgPCAxKSB7XG4gICAgcmV0dXJuIFtdXG4gIH1cblxuICAvL1NwZWNpYWwgY2FzZTogIEZvciAxRCB3ZSBjYW4ganVzdCBzb3J0IHRoZSBwb2ludHNcbiAgaWYoZCA9PT0gMSkge1xuICAgIHJldHVybiB0cmlhbmd1bGF0ZTFEKG4sIHBvaW50cywgaW5jbHVkZVBvaW50QXRJbmZpbml0eSlcbiAgfVxuICBcbiAgLy9MaWZ0IHBvaW50cywgc29ydFxuICB2YXIgbGlmdGVkID0gbmV3IEFycmF5KG4pXG4gIHZhciB1cHBlciA9IDEuMFxuICBmb3IodmFyIGk9MDsgaTxuOyArK2kpIHtcbiAgICB2YXIgcCA9IHBvaW50c1tpXVxuICAgIHZhciB4ID0gbmV3IEFycmF5KGQrMSlcbiAgICB2YXIgbCA9IDAuMFxuICAgIGZvcih2YXIgaj0wOyBqPGQ7ICsraikge1xuICAgICAgdmFyIHYgPSBwW2pdXG4gICAgICB4W2pdID0gdlxuICAgICAgbCArPSB2ICogdlxuICAgIH1cbiAgICB4W2RdID0gbFxuICAgIGxpZnRlZFtpXSA9IG5ldyBMaWZ0ZWRQb2ludCh4LCBpKVxuICAgIHVwcGVyID0gTWF0aC5tYXgobCwgdXBwZXIpXG4gIH1cbiAgdW5pcShsaWZ0ZWQsIGNvbXBhcmVMaWZ0ZWQpXG4gIFxuICAvL0RvdWJsZSBwb2ludHNcbiAgbiA9IGxpZnRlZC5sZW5ndGhcblxuICAvL0NyZWF0ZSBuZXcgbGlzdCBvZiBwb2ludHNcbiAgdmFyIGRwb2ludHMgPSBuZXcgQXJyYXkobiArIGQgKyAxKVxuICB2YXIgZGluZGV4ID0gbmV3IEFycmF5KG4gKyBkICsgMSlcblxuICAvL0FkZCBzdGVpbmVyIHBvaW50cyBhdCB0b3BcbiAgdmFyIHUgPSAoZCsxKSAqIChkKzEpICogdXBwZXJcbiAgdmFyIHkgPSBuZXcgQXJyYXkoZCsxKVxuICBmb3IodmFyIGk9MDsgaTw9ZDsgKytpKSB7XG4gICAgeVtpXSA9IDAuMFxuICB9XG4gIHlbZF0gPSB1XG5cbiAgZHBvaW50c1swXSA9IHkuc2xpY2UoKVxuICBkaW5kZXhbMF0gPSAtMVxuXG4gIGZvcih2YXIgaT0wOyBpPD1kOyArK2kpIHtcbiAgICB2YXIgeCA9IHkuc2xpY2UoKVxuICAgIHhbaV0gPSAxXG4gICAgZHBvaW50c1tpKzFdID0geFxuICAgIGRpbmRleFtpKzFdID0gLTFcbiAgfVxuXG4gIC8vQ29weSByZXN0IG9mIHRoZSBwb2ludHMgb3ZlclxuICBmb3IodmFyIGk9MDsgaTxuOyArK2kpIHtcbiAgICB2YXIgaCA9IGxpZnRlZFtpXVxuICAgIGRwb2ludHNbaSArIGQgKyAxXSA9IGgucG9pbnRcbiAgICBkaW5kZXhbaSArIGQgKyAxXSA9ICBoLmluZGV4XG4gIH1cblxuICAvL0NvbnN0cnVjdCBjb252ZXggaHVsbFxuICB2YXIgaHVsbCA9IGNoKGRwb2ludHMsIGZhbHNlKVxuICBpZihpbmNsdWRlUG9pbnRBdEluZmluaXR5KSB7XG4gICAgaHVsbCA9IGh1bGwuZmlsdGVyKGZ1bmN0aW9uKGNlbGwpIHtcbiAgICAgIHZhciBjb3VudCA9IDBcbiAgICAgIGZvcih2YXIgaj0wOyBqPD1kOyArK2opIHtcbiAgICAgICAgdmFyIHYgPSBkaW5kZXhbY2VsbFtqXV1cbiAgICAgICAgaWYodiA8IDApIHtcbiAgICAgICAgICBpZigrK2NvdW50ID49IDIpIHtcbiAgICAgICAgICAgIHJldHVybiBmYWxzZVxuICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgICAgICBjZWxsW2pdID0gdlxuICAgICAgfVxuICAgICAgcmV0dXJuIHRydWVcbiAgICB9KVxuICB9IGVsc2Uge1xuICAgIGh1bGwgPSBodWxsLmZpbHRlcihmdW5jdGlvbihjZWxsKSB7XG4gICAgICBmb3IodmFyIGk9MDsgaTw9ZDsgKytpKSB7XG4gICAgICAgIHZhciB2ID0gZGluZGV4W2NlbGxbaV1dXG4gICAgICAgIGlmKHYgPCAwKSB7XG4gICAgICAgICAgcmV0dXJuIGZhbHNlXG4gICAgICAgIH1cbiAgICAgICAgY2VsbFtpXSA9IHZcbiAgICAgIH1cbiAgICAgIHJldHVybiB0cnVlXG4gICAgfSlcbiAgfVxuXG4gIGlmKGQgJiAxKSB7XG4gICAgZm9yKHZhciBpPTA7IGk8aHVsbC5sZW5ndGg7ICsraSkge1xuICAgICAgdmFyIGggPSBodWxsW2ldXG4gICAgICB2YXIgeCA9IGhbMF1cbiAgICAgIGhbMF0gPSBoWzFdXG4gICAgICBoWzFdID0geFxuICAgIH1cbiAgfVxuXG4gIHJldHVybiBodWxsXG59IiwiXG5tb2R1bGUuZXhwb3J0cyA9IHBhcnNlXG5cbi8qKlxuICogZXhwZWN0ZWQgYXJndW1lbnQgbGVuZ3Roc1xuICogQHR5cGUge09iamVjdH1cbiAqL1xuXG52YXIgbGVuZ3RoID0ge2E6IDcsIGM6IDYsIGg6IDEsIGw6IDIsIG06IDIsIHE6IDQsIHM6IDQsIHQ6IDIsIHY6IDEsIHo6IDB9XG5cbi8qKlxuICogc2VnbWVudCBwYXR0ZXJuXG4gKiBAdHlwZSB7UmVnRXhwfVxuICovXG5cbnZhciBzZWdtZW50ID0gLyhbYXN0dnpxbWhsY10pKFteYXN0dnpxbWhsY10qKS9pZ1xuXG4vKipcbiAqIHBhcnNlIGFuIHN2ZyBwYXRoIGRhdGEgc3RyaW5nLiBHZW5lcmF0ZXMgYW4gQXJyYXlcbiAqIG9mIGNvbW1hbmRzIHdoZXJlIGVhY2ggY29tbWFuZCBpcyBhbiBBcnJheSBvZiB0aGVcbiAqIGZvcm0gYFtjb21tYW5kLCBhcmcxLCBhcmcyLCAuLi5dYFxuICpcbiAqIEBwYXJhbSB7U3RyaW5nfSBwYXRoXG4gKiBAcmV0dXJuIHtBcnJheX1cbiAqL1xuXG5mdW5jdGlvbiBwYXJzZShwYXRoKSB7XG5cdHZhciBkYXRhID0gW11cblx0cGF0aC5yZXBsYWNlKHNlZ21lbnQsIGZ1bmN0aW9uKF8sIGNvbW1hbmQsIGFyZ3Mpe1xuXHRcdHZhciB0eXBlID0gY29tbWFuZC50b0xvd2VyQ2FzZSgpXG5cdFx0YXJncyA9IHBhcnNlVmFsdWVzKGFyZ3MpXG5cblx0XHQvLyBvdmVybG9hZGVkIG1vdmVUb1xuXHRcdGlmICh0eXBlID09ICdtJyAmJiBhcmdzLmxlbmd0aCA+IDIpIHtcblx0XHRcdGRhdGEucHVzaChbY29tbWFuZF0uY29uY2F0KGFyZ3Muc3BsaWNlKDAsIDIpKSlcblx0XHRcdHR5cGUgPSAnbCdcblx0XHRcdGNvbW1hbmQgPSBjb21tYW5kID09ICdtJyA/ICdsJyA6ICdMJ1xuXHRcdH1cblxuXHRcdHdoaWxlICh0cnVlKSB7XG5cdFx0XHRpZiAoYXJncy5sZW5ndGggPT0gbGVuZ3RoW3R5cGVdKSB7XG5cdFx0XHRcdGFyZ3MudW5zaGlmdChjb21tYW5kKVxuXHRcdFx0XHRyZXR1cm4gZGF0YS5wdXNoKGFyZ3MpXG5cdFx0XHR9XG5cdFx0XHRpZiAoYXJncy5sZW5ndGggPCBsZW5ndGhbdHlwZV0pIHRocm93IG5ldyBFcnJvcignbWFsZm9ybWVkIHBhdGggZGF0YScpXG5cdFx0XHRkYXRhLnB1c2goW2NvbW1hbmRdLmNvbmNhdChhcmdzLnNwbGljZSgwLCBsZW5ndGhbdHlwZV0pKSlcblx0XHR9XG5cdH0pXG5cdHJldHVybiBkYXRhXG59XG5cbmZ1bmN0aW9uIHBhcnNlVmFsdWVzKGFyZ3Mpe1xuXHRhcmdzID0gYXJncy5tYXRjaCgvLT9bLjAtOV0rKD86ZVstK10/XFxkKyk/L2lnKVxuXHRyZXR1cm4gYXJncyA/IGFyZ3MubWFwKE51bWJlcikgOiBbXVxufVxuIiwiJ3VzZSBzdHJpY3QnO1xuXG52YXIgc2lnbiA9IHJlcXVpcmUoJy4vdXRpbGl0aWVzLmpzJykuc2lnbjtcbnZhciBjYWxjdWxhdGVEaXN0YW5jZSA9IHJlcXVpcmUoJy4vdXRpbGl0aWVzLmpzJykuZGlzdGFuY2U7XG5cbi8vIHZhciBwb2ludHMgPSByZXF1aXJlKCcuL2luaXRpYWxpemVQb2ludHMuanMnKS5wb2ludHM7XG4vLyB2YXIgY2l0eVNldCA9IHJlcXVpcmUoJy4vaW5pdGlhbGl6ZVBvaW50cy5qcycpLmNpdHlTZXQ7XG4vLyB2YXIgdGV4dFBvaW50c0lkID0gcmVxdWlyZSgnLi9pbml0aWFsaXplUG9pbnRzLmpzJykudGV4dFBvaW50c0lkO1xuLy8gdmFyIHBvc3NpYmxlU3RhcnRQb2ludHNJZCA9IHJlcXVpcmUoJy4vaW5pdGlhbGl6ZVBvaW50cy5qcycpLnBvc3NpYmxlU3RhcnRQb2ludHNJZDtcblxudmFyIGxpdmVNb3VzZVBvc2l0aW9uID0gcmVxdWlyZSgnLi9tb3VzZS5qcycpO1xuXG52YXIgVmVjdG9yID0gcmVxdWlyZSgnLi92ZWN0b3IuanMnKTtcblxudmFyIHJhbmRvbSA9IE1hdGgucmFuZG9tO1xudmFyIGZsb29yID0gTWF0aC5mbG9vcjtcbi8vIHZhciBSRVBVTFNJT04gPSAwLjA1O1xuLy8gdmFyIFJFUFVMU0lPTlNQRUVEID0gMC4wMDI7XG4vLyB2YXIgQU5UVkVMT0NJVFkgPSAwLjAwMTtcblxubW9kdWxlLmV4cG9ydHMgPSBmdW5jdGlvbihjb250YWluZXIsIGluaXRQb2ludHMsIG9wdGlvbnMpe1xuXG4gICAgY29uc29sZS5sb2coJ09wdGlvbnMgYW50IDonLCBvcHRpb25zKTtcbiAgICAvLyBEZWZpbmUgdGhvc2UgcGFyYW1ldGVycyBhcyBhdHRyaWJ1dGVzIG9mIEFudCBvYmplY3QgP1xuICAgIHZhciBSRVBVTFNJT04gPSBvcHRpb25zLnJlcFNpemU7XG4gICAgdmFyIFJFUFVMU0lPTlNQRUVEID0gb3B0aW9ucy5yZXBTcGVlZDtcbiAgICB2YXIgQU5UVkVMT0NJVFkgPSBvcHRpb25zLnZlbG9jaXR5O1xuICAgIHZhciBXRUlHSFQgPSBvcHRpb25zLndlaWdodDtcblxuICAgIHZhciBtb3VzZSA9IGxpdmVNb3VzZVBvc2l0aW9uKGNvbnRhaW5lcik7XG5cbiAgICB2YXIgcG9pbnRzID0gaW5pdFBvaW50cy5wb2ludHM7XG4gICAgdmFyIGNpdHlTZXQgPSBpbml0UG9pbnRzLmNpdHlTZXQ7XG4gICAgdmFyIHRleHRQb2ludHNJZCA9IGluaXRQb2ludHMudGV4dFBvaW50c0lkO1xuICAgIHZhciBwb3NzaWJsZVN0YXJ0UG9pbnRzSWQgPSBpbml0UG9pbnRzLnBvc3NpYmxlU3RhcnRQb2ludHNJZDtcblxuXG4gICAgZnVuY3Rpb24gQW50KHBvaW50KSB7XG4gICAgICAgIHRoaXMueCA9IHBvaW50Lng7ICAgICAgICAgICAgICAgIFxuICAgICAgICB0aGlzLnkgPSBwb2ludC55O1xuICAgICAgICB0aGlzLnZlbG9jaXR5ID0gQU5UVkVMT0NJVFk7XG4gICAgICAgIHRoaXMud2VpZ2h0ID0gV0VJR0hUO1xuICAgICAgICB0aGlzLnJlcFNpemUgPSBSRVBVTFNJT047XG4gICAgICAgIHRoaXMucmVwU3BlZWQgPSBSRVBVTFNJT05TUEVFRDtcbiAgICAgICAgdGhpcy5lZGdlID0gdW5kZWZpbmVkO1xuICAgICAgICB0aGlzLnN0YXRlID0gXCJmb3JhZ2VcIjtcbiAgICAgICAgdGhpcy5lZGdlcyA9IFtdO1xuICAgICAgICB0aGlzLmxhc3RDaXR5ID0gdW5kZWZpbmVkO1xuICAgICAgICB0aGlzLm9yaWdpbiA9IHBvaW50O1xuICAgICAgICB0aGlzLmRlc3RpbmF0aW9uID0gdW5kZWZpbmVkO1xuICAgICAgICB0aGlzLm9yaWVudGF0aW9uID0gdW5kZWZpbmVkO1xuICAgICAgICB0aGlzLmRpcmVjdGlvbiA9IG5ldyBWZWN0b3IoMCwwKTtcbiAgICAgICAgdGhpcy5wcm9nID0gMDtcbiAgICB9XG4gICAgLy8gZm9yYWdlOiB0aGUgYW50IHdhbmRlcnMgYXJvdW5kIHdpdGhvdXQgYW55IHBoZXJvbW9uIGRlcG9zaXRpb25cbiAgICAvLyBvbmNlIGl0IGZpbmRzIGEgY2l0eSwgaXQgc3RhcnRzIHJlbWVtYmVyaW5nIHRoZSBub2RlcyBpdCBnb2VzIHRocm91Z2hcbiAgICAvLyB3aGVuIGl0IGZpbmRzIGFub3RoZXIgY2l0eSwgaXQgY29tcHV0ZXMgdGhlIHBhdGggbGVuZ3RoIGFuZCBhZGRzIHBoZXJvbW9ucyBvbmUgZWFjaCBlZGdlc1xuICAgIC8vIHByb3BvcnRpb25uYWx5IHRvIHRoZSBzaG9ydGVzdG5lc3Mgb2YgdGhlIHBhdGhcbiAgICAvLyBpdCByZXNldHMgdGhlIGxpc3Qgb2Ygbm9kZXMgYW5kIGNvbnRpbnVlc1xuICAgIC8vIHdoaWxlIGZvcmFnaW5nIHRoZSBhbnQgY2hvc2VzIHRoZSBwYXRoIHdpdGggYSBwaGVyb21vbiBwcmVmZXJlbmNlXG5cblxuICAgIC8vIHN0YXRpYyBtZXRob2RzXG4gICAgQW50LmdlbmVyYXRlUmFuZFN0YXJ0UG9pbnQgPSBmdW5jdGlvbigpIHtcbiAgICAgICAgdmFyIHJhbmRJZCA9IE1hdGguZmxvb3IocG9zc2libGVTdGFydFBvaW50c0lkLmxlbmd0aCAqIHJhbmRvbSgpKTtcbiAgICAgICAgdmFyIHJhbmRTdGFydFBvaW50ID0gcG9pbnRzW3Bvc3NpYmxlU3RhcnRQb2ludHNJZFtyYW5kSWRdXTtcbiAgICAgICAgcmV0dXJuIHJhbmRTdGFydFBvaW50O1xuICAgIH1cblxuXG4gICAgLy8gbWV0aG9kc1xuICAgIEFudC5wcm90b3R5cGUgPSB7XG5cbiAgICAgICAgdHJhbnNpdDogZnVuY3Rpb24oKXtcbiAgICAgICAgICAgIHN3aXRjaCAodGhpcy5zdGF0ZSkge1xuICAgICAgICAgICAgY2FzZSBcImZvcmFnZVwiOlxuICAgICAgICAgICAgICAgIHZhciByZXMgPSB0aGlzLm1vdmUoKTtcbiAgICAgICAgICAgICAgICBpZiAocmVzLmNpdHlSZWFjaGVkKSB7XG4gICAgICAgICAgICAgICAgICAgIHRoaXMuc3RhdGUgPSBcInBoZXJvbW9uaW5nXCI7XG4gICAgICAgICAgICAgICAgICAgIHRoaXMubGFzdENpdHkgPSB0aGlzLm9yaWdpbi5pZDtcbiAgICAgICAgICAgICAgICB9O1xuICAgICAgICAgICAgICAgIGJyZWFrO1xuICAgICAgICAgICAgY2FzZSBcInBoZXJvbW9uaW5nXCI6XG4gICAgICAgICAgICAgICAgdmFyIHJlcyA9IHRoaXMubW92ZSgpO1xuICAgICAgICAgICAgICAgIGlmIChyZXMuZWRnZUNoYW5nZWQpIHtcbiAgICAgICAgICAgICAgICAgICAgdGhpcy5lZGdlcy5wdXNoKHRoaXMuZWRnZSk7XG4gICAgICAgICAgICAgICAgICAgIC8vIGZvdW5kIGEgY2l0eVxuICAgICAgICAgICAgICAgICAgICBpZiAocmVzLmNpdHlSZWFjaGVkICYmICh0aGlzLm9yaWdpbi5pZCAhPSB0aGlzLmxhc3RDaXR5KSApe1xuICAgICAgICAgICAgICAgICAgICAgICAgLy8gY29tcHV0ZSB0aGUgbGVuZ3RoIG9mIHRoZSBwYXRoXG4gICAgICAgICAgICAgICAgICAgICAgICB2YXIgcGF0aExlbmd0aCA9IHRoaXMuZWRnZXMubWFwKGZ1bmN0aW9uKGUpe3JldHVybiBlLmRpc3RhbmNlfSkucmVkdWNlKGZ1bmN0aW9uKGEsYil7cmV0dXJuIGEgKyBifSk7XG4gICAgICAgICAgICAgICAgICAgICAgICB2YXIgZGVsdGFQaGVyb21vbmUgPSAxL3BhdGhMZW5ndGg7XG4gICAgICAgICAgICAgICAgICAgICAgICB2YXIgYW50V2VpZ2h0ID0gdGhpcy53ZWlnaHQ7XG4gICAgICAgICAgICAgICAgICAgICAgICB0aGlzLmVkZ2VzLmZvckVhY2goZnVuY3Rpb24oZSl7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgdmFyIGEgPSBlLnB0MSwgYiA9IGUucHQyLCB3ZWlnaHQgPSAxOyAgXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgLy8gaW5jcmVhc2VkIGRyb3BwZWQgcGhlcm9tb25zIGZvciB0ZXh0RWRnZXNcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBpZiAoKGNpdHlTZXQuaW5kZXhPZihhLmlkKSAhPSAtMSkgJiYgY2l0eVNldC5pbmRleE9mKGIuaWQpICE9IC0xICYmIChNYXRoLmFicyhhLmlkIC0gYi5pZCkgPT0gMSkpXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB3ZWlnaHQgKj0gYW50V2VpZ2h0O1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBlLnBoZXJvbW9uICs9IChkZWx0YVBoZXJvbW9uZSAqIHdlaWdodCk7XG4gICAgICAgICAgICAgICAgICAgICAgICB9KTtcblxuICAgICAgICAgICAgICAgICAgICAgICAgdGhpcy5lZGdlcyA9IFt0aGlzLmVkZ2VdO1xuICAgICAgICAgICAgICAgICAgICAgICAgdGhpcy5sYXN0Q2l0eSA9IHRoaXMub3JpZ2luLmlkO1xuICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICBicmVhaztcbiAgICAgICAgICAgIH1cblxuICAgICAgICB9LFxuXG4gICAgICAgIHNldERpcmVjdGlvbjogZnVuY3Rpb24oKXtcbiAgICAgICAgICAgIHZhciBwb3NzaWJsZUVkZ2VzID0gW107XG5cbiAgICAgICAgICAgIGZvciAodmFyIGkgPSAwOyBpIDwgdGhpcy5vcmlnaW4ubmV4dHMubGVuZ3RoOyBpKyspXG4gICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgcG9zc2libGVFZGdlc1tpXSA9IHRoaXMub3JpZ2luLm5leHRzW2ldO1xuICAgICAgICAgICAgfSBcblxuICAgICAgICAgICAgLy8gY29uc29sZS5sb2coJ3NtZWxsczE6ICcsIHBvc3NpYmxlRWRnZXMpO1xuXG4gICAgICAgICAgICBwb3NzaWJsZUVkZ2VzLnNwbGljZShwb3NzaWJsZUVkZ2VzLmluZGV4T2YodGhpcy5lZGdlKSwxKTtcblxuICAgICAgICAgICAgLy8gZmxpcCBhIGNvaW4gYW5kIGVpdGhlciB0YWtlIHRoZSBzbWVsbGllc3QgcGF0aCBvciBhIHJhbmRvbSBvbmVcbiAgICAgICAgICAgIGlmIChyYW5kb20oKSA8IDAuOCl7XG4gICAgICAgICAgICAgICAgdmFyIHNtZWxscyA9IHBvc3NpYmxlRWRnZXMubWFwKGZ1bmN0aW9uKGUpe3JldHVybiBlLnBoZXJvbW9uO30pO1xuICAgICAgICAgICAgICAgIHZhciBpbmRleCA9IHNtZWxscy5pbmRleE9mKE1hdGgubWF4LmFwcGx5KE1hdGgsIHNtZWxscykpO1xuICAgICAgICAgICAgICAgIHRoaXMuZWRnZSA9IHBvc3NpYmxlRWRnZXNbaW5kZXhdO1xuICAgICAgICAgICAgfSBcbiAgICAgICAgICAgIGVsc2V7XG4gICAgICAgICAgICAgICAgdGhpcy5lZGdlID0gcG9zc2libGVFZGdlc1tmbG9vcihyYW5kb20oKSpwb3NzaWJsZUVkZ2VzLmxlbmd0aCldO1xuICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIFxuXG4gICAgICAgICAgICAvLyBzZXQgdGhlIGRlc3RpbmF0aW9uIHBvaW50LCBiZWluZyBlZGdlLnB0MSBvciBlZGdlLnB0MlxuICAgICAgICAgICAgdGhpcy5kZXN0aW5hdGlvbiA9ICh0aGlzLm9yaWdpbiA9PSB0aGlzLmVkZ2UucHQxKSA/IHRoaXMuZWRnZS5wdDIgOiB0aGlzLmVkZ2UucHQxO1xuXG4gICAgICAgICAgICB0aGlzLmRpcmVjdGlvbi54ID0gdGhpcy5kZXN0aW5hdGlvbi54IC0gdGhpcy5vcmlnaW4ueDsgXG4gICAgICAgICAgICB0aGlzLmRpcmVjdGlvbi55ID0gdGhpcy5kZXN0aW5hdGlvbi55IC0gdGhpcy5vcmlnaW4ueTtcblxuICAgICAgICAgICAgdGhpcy5kaXJlY3Rpb24ubm9ybWFsaXplKCk7XG4gICAgICAgIH0sXG5cbiAgICAgICAgbW92ZTogZnVuY3Rpb24oKXtcbiAgICAgICAgICAgIC8vIGNvbnNvbGUubG9nKCdtb3ZlJyk7XG4gICAgICAgICAgICB2YXIgZWRnZUNoYW5nZWQ7XG4gICAgICAgICAgICB2YXIgY2l0eVJlYWNoZWQgPSBmYWxzZTtcblxuICAgICAgICAgICAgdGhpcy5kaXJlY3Rpb24ueCA9IHRoaXMuZGVzdGluYXRpb24ueCAtIHRoaXMueDsgXG4gICAgICAgICAgICB0aGlzLmRpcmVjdGlvbi55ID0gdGhpcy5kZXN0aW5hdGlvbi55IC0gdGhpcy55O1xuICAgICAgICAgICAgdGhpcy5kaXJlY3Rpb24ubm9ybWFsaXplKCk7XG5cbiAgICAgICAgICAgIC8vIG9uIGVkZ2VcbiAgICAgICAgICAgIGlmICgoY2FsY3VsYXRlRGlzdGFuY2UodGhpcywgdGhpcy5kZXN0aW5hdGlvbikgPiB0aGlzLnJlcFNwZWVkKSl7XG5cbiAgICAgICAgICAgICAgICAvLyBhIGRlbHRhIG1vdmVtZW50IHdpbGwgYmUgYXBwbGllZCBpZiBjb2xsaXNpb24gd2l0aCBvYnN0YWNsZSBkZXRlY3RlZFxuICAgICAgICAgICAgICAgIHZhciBkZWx0YSA9IHRoaXMuYXZvaWRPYnN0YWNsZSgpO1xuXG4gICAgICAgICAgICAgICAgdGhpcy54ICs9IHRoaXMudmVsb2NpdHkgKiB0aGlzLmRpcmVjdGlvbi54ICsgZGVsdGEueCAqIHRoaXMucmVwU3BlZWQ7XG4gICAgICAgICAgICAgICAgdGhpcy55ICs9IHRoaXMudmVsb2NpdHkgKiB0aGlzLmRpcmVjdGlvbi55ICsgZGVsdGEueSAqIHRoaXMucmVwU3BlZWQ7XG5cbiAgICAgICAgICAgICAgICB0aGlzLnByb2cgPSB0aGlzLmNhbGN1bGF0ZVByb2dyZXNzaW9uKCk7XG4gICAgICAgICAgICAgICAgXG4gICAgICAgICAgICAgICAgZWRnZUNoYW5nZWQgPSBmYWxzZTtcblxuICAgICAgICAgICAgLy8gb24gdmVydGV4XG4gICAgICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgICAgIC8vIGNvbnNvbGUubG9nKCdyZWFjaGVkJyk7XG4gICAgICAgICAgICAgICAgdGhpcy5zdGVwID0gMDtcbiAgICAgICAgICAgICAgICB0aGlzLnByb2cgPSAwO1xuICAgICAgICAgICAgICAgIHRoaXMub3JpZ2luID0gdGhpcy5kZXN0aW5hdGlvbjtcbiAgICAgICAgICAgICAgICB0aGlzLnggPSB0aGlzLm9yaWdpbi54O1xuICAgICAgICAgICAgICAgIHRoaXMueSA9IHRoaXMub3JpZ2luLnk7XG5cbiAgICAgICAgICAgICAgICB0aGlzLnNldERpcmVjdGlvbigpO1xuXG4gICAgICAgICAgICAgICAgY2l0eVJlYWNoZWQgPSAoY2l0eVNldC5pbmRleE9mKHRoaXMub3JpZ2luLmlkKSAhPSAtMSk7XG4gICAgICAgICAgICAgICAgZWRnZUNoYW5nZWQgPSB0cnVlO1xuICAgICAgICAgICAgfVxuICAgICAgICAgICAgcmV0dXJuIHtjaXR5UmVhY2hlZDogY2l0eVJlYWNoZWQsIGVkZ2VDaGFuZ2VkOiBlZGdlQ2hhbmdlZH07XG4gICAgICAgIH0sXG5cbiAgICAgICAgYXZvaWRPYnN0YWNsZTogZnVuY3Rpb24oKXtcbiAgICAgICAgICAgIHZhciBkaXN0YW5jZSA9IGNhbGN1bGF0ZURpc3RhbmNlKHRoaXMsIG1vdXNlKTtcbiAgICAgICAgXG4gICAgICAgICAgICBpZiAoZGlzdGFuY2UgPD0gdGhpcy5yZXBTaXplKSB7XG5cbiAgICAgICAgICAgICAgICByZXR1cm4ge1xuICAgICAgICAgICAgICAgICAgICAvLyBkZWx0YSBtb3ZlbWVudCBpcyBjb21wb3NlZCBvZiBhIHJlcHVsc2lvbiBkZWx0YSBhbmQgYSBjaXJjdWxhciBkZWx0YSBcbiAgICAgICAgICAgICAgICAgICAgeDogKHRoaXMueCAtIG1vdXNlLngpL2Rpc3RhbmNlICsgKHRoaXMueSAtIG1vdXNlLnkpL2Rpc3RhbmNlICogMSxcbiAgICAgICAgICAgICAgICAgICAgeTogKHRoaXMueSAtIG1vdXNlLnkpL2Rpc3RhbmNlIC0gKHRoaXMueCAtIG1vdXNlLngpL2Rpc3RhbmNlICogMVxuICAgICAgICAgICAgICAgIH07XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICBlbHNlXG4gICAgICAgICAgICAgICAgcmV0dXJuIHt4OjAsIHk6MH07XG4gICAgICAgIH0sXG5cbiAgICAgICAgY2FsY3VsYXRlUHJvZ3Jlc3Npb246IGZ1bmN0aW9uKCl7XG4gICAgICAgICAgICB2YXIgdiA9IG5ldyBWZWN0b3IodGhpcy54IC0gdGhpcy5vcmlnaW4ueCwgdGhpcy55IC0gdGhpcy5vcmlnaW4ueSk7XG4gICAgICAgICAgICB2YXIgbm9ybSA9IHYubm9ybSgpO1xuXG4gICAgICAgICAgICB2YXIgdGhldGEgPSAodi54ICogdGhpcy5lZGdlLmRpcmVjdGlvbi54ICsgdi55ICogdGhpcy5lZGdlLmRpcmVjdGlvbi55KSAvIG5vcm07XG4gICAgICAgICAgICB2YXIgcHJvZyA9IG5vcm0gKiBNYXRoLmFicyh0aGV0YSk7XG4gICAgICAgICAgICAvLyByZXR1cm5zIGxlbmd0aCBvZiBwcm9qZWN0aW9uIG9uIGVkZ2VcbiAgICAgICAgICAgIHJldHVybiBwcm9nO1xuICAgICAgICB9XG5cbiAgICB9O1xuICAgIHJldHVybiBBbnQ7XG59XG5cbiIsIid1c2Ugc3RyaWN0J1xuXG5tb2R1bGUuZXhwb3J0cyA9IGZ1bmN0aW9uIChBbnQpIHtcblxuXHR2YXIgbmJBbnRzUGVyU3RlcCA9IDEwMDtcblxuXHRmdW5jdGlvbiBjcmVhdGVHcm91cChwb3B1bGF0aW9uKXtcblx0XHRmb3IgKHZhciBpID0gMDsgaSA8IG5iQW50c1BlclN0ZXA7IGkrKykge1xuXHRcdFx0dmFyIG5ld0FudCA9IG5ldyBBbnQoQW50LmdlbmVyYXRlUmFuZFN0YXJ0UG9pbnQoKSk7XG5cdFx0XHRuZXdBbnQuc2V0RGlyZWN0aW9uKCk7XG5cdFx0XHRwb3B1bGF0aW9uLnB1c2gobmV3QW50KTtcblx0XHR9XG5cbi8vIFx0XHRjb25zb2xlLmxvZygnQ3JlYXRlZCBBbnRzIEdyb3VwOiBcXFxuLy8gKCsgJyArIG5iQW50c1BlclN0ZXAgKyAnKSA9PiAnICsgcG9wdWxhdGlvbi5sZW5ndGgpO1xuXG5cdFx0cmV0dXJuIHBvcHVsYXRpb247XG5cdH1cblxuXHRmdW5jdGlvbiByZW1vdmVHcm91cChwb3B1bGF0aW9uLCBuYkRlYWQpe1xuXHRcdHBvcHVsYXRpb24gPSBwb3B1bGF0aW9uLnNsaWNlKDAsIHBvcHVsYXRpb24ubGVuZ3RoIC0gbmJEZWFkKTtcblxuLy8gXHRcdGNvbnNvbGUubG9nKCdSZW1vdmVkIEFudHMgR3JvdXA6IFxcXG4vLyAoLSAnICsgbmJBbnRzUGVyU3RlcCArICcpID0+ICcgKyBwb3B1bGF0aW9uLmxlbmd0aCk7XG5cblx0XHRyZXR1cm4gcG9wdWxhdGlvbjtcblxuXHR9XG5cblx0cmV0dXJuIHtcblx0XHRjcmVhdGU6IGNyZWF0ZUdyb3VwLFxuXHRcdHJlbW92ZTogcmVtb3ZlR3JvdXBcblx0fTtcblxufVxuXHQiLCIndXNlIHN0cmljdCdcblxudmFyIGR0ID0gcmVxdWlyZShcImRlbGF1bmF5LXRyaWFuZ3VsYXRlXCIpO1xuXG52YXIgcmFuZ2UgPSByZXF1aXJlKCcuL3V0aWxpdGllcy5qcycpLnJhbmdlO1xuXG4vLyB2YXIgcG9pbnRzID0gcmVxdWlyZSgnLi9pbml0aWFsaXplUG9pbnRzLmpzJykucG9pbnRzO1xudmFyIHRleHRNZXNoID0gcmVxdWlyZSgnLi9pbml0aWFsaXplUG9pbnRzLmpzJykudGV4dE1lc2g7XG52YXIgY2l0eVNldCA9IHJlcXVpcmUoJy4vaW5pdGlhbGl6ZVBvaW50cy5qcycpLmNpdHlTZXQ7XG52YXIgbmJSYW5kb21Qb2ludHMgPSByZXF1aXJlKCcuL2luaXRpYWxpemVQb2ludHMuanMnKS5uYlJhbmRvbVBvaW50cztcbnZhciBmb3JjZWRFZGdlcyA9IHJlcXVpcmUoJy4vaW5pdGlhbGl6ZVBvaW50cy5qcycpLmZvcmNlZEVkZ2VzO1xuXG52YXIgRWRnZSA9IHJlcXVpcmUoJy4vZWRnZS5qcycpO1xuXG5cbm1vZHVsZS5leHBvcnRzID0gZnVuY3Rpb24ocG9pbnRzKXtcbiAgICAvLyB0cmlhbmd1bGF0ZVxuICAgIHZhciBjZWxscyA9IGR0KHBvaW50cy5tYXAoZnVuY3Rpb24ocCl7XG4gICAgICAgIHJldHVybiBbcC54LCBwLnldO1xuICAgIH0pKTtcblxuICAgIHZhciBlZGdlcyA9IFtdO1xuICAgIHZhciBwZXJtdXRhdGlvbnMgPSBbWzAsMV0sIFswLDJdLCBbMSwyXV07XG5cbiAgICAvLyBmb3JjZSB0aGUgZWRnZXMgb2YgdGhlIHRleHQgdG8gYmUgZWRnZXMgb2YgdGhlIGdyYXBoXG4gICAgaWYgKHRleHRNZXNoKSB7XG4gICAgICAgIHJhbmdlKDAsIHBvaW50cy5sZW5ndGggLSBuYlJhbmRvbVBvaW50cykuZm9yRWFjaChmdW5jdGlvbihpZCl7XG4gICAgICAgICAgICB2YXIgZGlyZWN0TGluayA9IGZvcmNlZEVkZ2VzW2lkXTtcbiAgICAgICAgICAgIHZhciB0ZXh0RWRnZSA9IEVkZ2UuY3JlYXRlKHBvaW50c1tpZF0sIHBvaW50c1tkaXJlY3RMaW5rXSk7XG4gICAgICAgICAgICBlZGdlcy5wdXNoKHRleHRFZGdlKTtcbiAgICAgICAgICAgIHBvaW50c1tpZF0ubmV4dHMucHVzaCh0ZXh0RWRnZSk7XG4gICAgICAgIH0pXG4gICAgfVxuXG5cbiAgICBjZWxscy5mb3JFYWNoKGZ1bmN0aW9uKGNlbGwpe1xuICAgICAgIFxuICAgICAgICBmb3IgKHZhciBpID0gMDsgaSA8IDM7ICsraSl7ICAvLyBmb3IgZWFjaCBwb2ludC5pZCBsaXN0ZWQgaW4gY3VycmVudCBjZWxsXG4gICAgICAgICAgICB2YXIgcHQgPSBwb2ludHNbY2VsbFtpXV07XG5cbiAgICAgICAgICAgIGZvciAodmFyIGogPSAxOyBqIDwgMzsgKytqKXsgXG5cbiAgICAgICAgICAgICAgICB2YXIgcHRqID0gcG9pbnRzW2NlbGxbKCBpICsgaiApICUgM11dOyAvLyBwaWNrIG9uZSBvZiB0aGUgb3RoZXIgMiBwb2ludHMgb2YgdGhlIGNlbGxcbiAgICAgICAgICAgICAgICB2YXIgbmV3RWRnZSA9IHVuZGVmaW5lZDtcblxuICAgICAgICAgICAgICAgIC8vIGlmIHB0IGFscmVhZHkgaGFzIG5leHRFZGdlcyAuLi5cbiAgICAgICAgICAgICAgICBpZiAocHQubmV4dHMubGVuZ3RoICE9IDApIHtcbiAgICAgICAgICAgICAgICAgICAgXG4gICAgICAgICAgICAgICAgICAgIC8vIC4uLiBnZXQgdGhlIHBvaW50cyBjb3JyZXNwb25kaW5nIC4uLlxuICAgICAgICAgICAgICAgICAgICB2YXIgdGVtcFBvaW50cyA9IHB0Lm5leHRzLm1hcChmdW5jdGlvbihlKXtcbiAgICAgICAgICAgICAgICAgICAgICAgIHJldHVybiBbZS5wdDEsIGUucHQyXTtcbiAgICAgICAgICAgICAgICAgICAgfSkucmVkdWNlKGZ1bmN0aW9uKGEsIGIpe1xuICAgICAgICAgICAgICAgICAgICAgICAgIHJldHVybiBhLmNvbmNhdChiKTtcbiAgICAgICAgICAgICAgICAgICAgfSk7XG5cbiAgICAgICAgICAgICAgICAgICAgLy8gLi4uIGFuZCBjaGVjayBpZiBwdGogYWxyZWFkeSBpcyBwYXJ0IG9mIHRoZSBleGlzdGluZyBuZXh0RWRnZXMuIElmIG5vdCwgYWRkIHRoZSBlZGdlLlxuICAgICAgICAgICAgICAgICAgICBpZiAodGVtcFBvaW50cy5pbmRleE9mKHB0aikgPT0gLTEpe1xuICAgICAgICAgICAgICAgICAgICAgICAgbmV3RWRnZSA9IEVkZ2UuY3JlYXRlKHB0LCBwdGopO1xuICAgICAgICAgICAgICAgICAgICAgICAgZWRnZXMucHVzaChuZXdFZGdlKTtcbiAgICAgICAgICAgICAgICAgICAgICAgIHB0Lm5leHRzLnB1c2gobmV3RWRnZSk7XG4gICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgZWxzZSB7XG4gICAgICAgICAgICAgICAgICAgIG5ld0VkZ2UgPSBFZGdlLmNyZWF0ZShwdCwgcHRqKTtcbiAgICAgICAgICAgICAgICAgICAgZWRnZXMucHVzaChuZXdFZGdlKTtcbiAgICAgICAgICAgICAgICAgICAgcHQubmV4dHMucHVzaChuZXdFZGdlKTtcbiAgICAgICAgICAgICAgICB9XG5cbiAgICAgICAgICAgICAgICAvLyBhZGQgYWxzbyB0aGUgZWRnZSB0byB0aGUgZWRnZSdzIG90aGVyIHBvaW50J3MgbmV4dEVkZ2VzXG4gICAgICAgICAgICAgICAgaWYgKG5ld0VkZ2UgIT0gdW5kZWZpbmVkKXtcbiAgICAgICAgICAgICAgICAgICAgcHRqLm5leHRzLnB1c2gobmV3RWRnZSk7XG4gICAgICAgICAgICAgICAgfSAgICAgICAgIFxuICAgICAgICAgICAgfVxuXG4gICAgICAgICAgICAvLyBhZGQgdGhlIHRleHRFZGdlcyB0byBuZXh0RWRnZXMgbWFwXG4gICAgICAgICAgICBpZiAodGV4dE1lc2ggJiYgKGNpdHlTZXQuaW5kZXhPZihwdCkgIT0gLTEpKSB7XG4gICAgICAgICAgICAgICAgdmFyIHRleHRFZGdlID0gRWRnZS5jcmVhdGUocHQsIHBvaW50c1twdC5pZCArIDFdKTtcbiAgICAgICAgICAgICAgICBlZGdlcy5wdXNoKHRleHRFZGdlKTtcbiAgICAgICAgICAgICAgICBwdC5uZXh0cy5wdXNoKHRleHRFZGdlKTtcbiAgICAgICAgICAgIH1cblxuICAgICAgICB9XG4gICAgfSk7XG5cbiAgICByZXR1cm4gZWRnZXM7XG59OyIsIid1c2Ugc3RyaWN0JztcblxudmFyIHNxcnQgPSBNYXRoLnNxcnQ7XG52YXIgcG93ID0gTWF0aC5wb3c7XG52YXIgYWJzID0gTWF0aC5hYnM7XG52YXIgYXRhbiA9IE1hdGguYXRhbjtcblxudmFyIFZlY3RvciA9IHJlcXVpcmUoJy4vdmVjdG9yLmpzJyk7XG5cblxuZnVuY3Rpb24gRWRnZShwdEEsIHB0Qikge1xuICAgIHZhciBkaXN0YW5jZSA9IHNxcnQoIHBvdyhwdEEueCAtIHB0Qi54LCAyKSArIHBvdyhwdEEueSAtIHB0Qi55LCAyKSApO1xuXG4gICAgLy8gZmluZCBsaW5lIGVxdWF0aW9uIGF4ICsgYnkgKyBjID0gMFxuICAgIHZhciBhID0gMTtcbiAgICB2YXIgYiA9IC0gKHB0Qi54IC0gcHRBLngpIC8gKHB0Qi55IC0gcHRBLnkpO1xuXG4gICAgLy8gb3JpZW50YXRlIHZlY3RvciAoYSxiKVxuICAgIGlmIChiIDwgMCl7XG4gICAgICAgIGIgPSAtYjtcbiAgICAgICAgYSA9IC1hO1xuICAgIH1cblxuICAgIC8vIG5vcm1hbGl6ZSB2ZWN0b3IgKGEsYilcbiAgICB2YXIgbiA9IG5ldyBWZWN0b3IoYSwgYik7XG4gICAgbi5ub3JtYWxpemUoKTtcblxuICAgIHZhciBjID0gLSAoYSAqIHB0QS54ICsgYiAqIHB0QS55KTtcblxuICAgIC8vIC8vIGNhbGN1bGF0ZSB2ZWN0b3IgZGlyZWN0b3JcbiAgICB2YXIgdiA9IG5ldyBWZWN0b3IocHRCLnggLSBwdEEueCwgcHRCLnkgLSBwdEEueSk7XG4gICAgXG4gICAgdi5ub3JtYWxpemUoKTtcblxuICAgIHRoaXMuaWQgPSB1bmRlZmluZWQ7XG4gICAgdGhpcy5wdDEgPSBwdEE7XG4gICAgdGhpcy5wdDIgPSBwdEI7XG4gICAgdGhpcy5kaXJlY3Rpb24gPSB2O1xuICAgIHRoaXMub3J0aERpcmVjdGlvbiA9IG47IFxuICAgIHRoaXMuZGlzdGFuY2UgPSBkaXN0YW5jZTtcbiAgICB0aGlzLnBoZXJvbW9uID0gMS9kaXN0YW5jZTtcbiAgICB0aGlzLmxpbmUgPSB7XG4gICAgICAgIGE6IGEsXG4gICAgICAgIGI6IGIsXG4gICAgICAgIGM6IGMsXG4gICAgfTtcblxuICAgIGlmICh0aGlzLmRpc3RhbmNlID09PSAwKSBjb25zb2xlLmxvZygnWkVSTyAhJyk7XG59XG5cblxuLy8gc3RhdGljIG1ldGhvZHNcbkVkZ2UuY3JlYXRlID0gZnVuY3Rpb24ocHRBLCBwdEIpIHtcbiAgICB2YXIgZWRnZSA9IG5ldyBFZGdlKHB0QSwgcHRCKTtcbiAgICByZXR1cm4gZWRnZTtcbn1cblxuXG4vLyBtZXRob2RzXG5FZGdlLnByb3RvdHlwZSA9IHtcblxuICAgIGdldE90aGVyUG9pbnQ6IGZ1bmN0aW9uKHBvaW50KSB7XG4gICAgICAgIGlmIChwb2ludCA9PSB0aGlzLnB0MSlcbiAgICAgICAgICAgIHJldHVybiB0aGlzLnB0MjtcbiAgICAgICAgZWxzZSBpZiAocG9pbnQgPT0gdGhpcy5wdDIpXG4gICAgICAgICAgICByZXR1cm4gdGhpcy5wdDE7XG4gICAgICAgIGVsc2VcbiAgICAgICAgICAgIGNvbnNvbGUubG9nKFwiRXJyb3JcIik7XG4gICAgfSxcblxuICAgIGNhbGN1bGF0ZURpc3RhbmNlOiBmdW5jdGlvbih4LCB5KSB7XG4gICAgICAgIHZhciBhID0gdGhpcy5saW5lLmEsXG4gICAgICAgICAgICBiID0gdGhpcy5saW5lLmIsXG4gICAgICAgICAgICBjID0gdGhpcy5saW5lLmM7XG4gICAgICAgIHJldHVybiBhYnMoYSAqIHggKyBiICogeSArIGMpIC8gTWF0aC5zcXJ0KE1hdGgucG93KGEsMikgKyBNYXRoLnBvdyhiLDIpKTtcbiAgICB9LFxuXG59XG5tb2R1bGUuZXhwb3J0cyA9IEVkZ2U7IiwiJ3VzZSBzdHJpY3QnO1xuXG52YXIgcGFyc2UgPSByZXF1aXJlKCdwYXJzZS1zdmctcGF0aCcpO1xuXG52YXIgcmFuZ2UgPSByZXF1aXJlKCcuL3V0aWxpdGllcy5qcycpLnJhbmdlO1xuXG52YXIgUG9pbnQgPSByZXF1aXJlKCcuL3BvaW50LmpzJyk7XG52YXIgc3ZnUGF0aCA9IHJlcXVpcmUoJy4vc3ZnUGF0aC5qcycpO1xuXG52YXIgcmFuZG9tID0gTWF0aC5yYW5kb207XG5cbnZhciBuYkNpdHkgPSAyO1xuXG52YXIgdGV4dE1lc2ggPSB0cnVlO1xuXG4vLyBGcmFtZSBkZWZpbml0aW9uXG52YXIgeEluaXQgPSAwLCB5SW5pdCA9IDA7XG52YXIgdyA9IDEsXG4gICAgaCA9IDE7XG5cbnZhciBzdmdTdHJpbmcgPSBzdmdQYXRoO1xuXG5mdW5jdGlvbiBzdmdUb1BvaW50cyhzdmdTdHJpbmcpIHtcbiAgICB2YXIgcG9pbnRzID0gW107XG4gICAgdmFyIGVkZ2VzID0gT2JqZWN0LmNyZWF0ZShudWxsKTtcblxuICAgIHZhciBiZWdpbmluZ1BhdGg7XG5cbiAgICB2YXIgWCA9IDA7XG4gICAgdmFyIFkgPSAwO1xuICAgIHZhciBuYlBvaW50cyA9IDA7XG4gICAgdmFyIHByZXZQb2ludDtcblxuICAgIHZhciBjb21tYW5kcyA9IHBhcnNlKHN2Z1N0cmluZylcbiAgICBmb3IgKHZhciBpPTA7IGk8Y29tbWFuZHMubGVuZ3RoOyBpKyspe1xuICAgICAgICB2YXIgY29tbWFuZCA9IGNvbW1hbmRzW2ldO1xuICAgICAgICBzd2l0Y2ggKGNvbW1hbmRbMF0pIHtcbiAgICAgICAgICAgIGNhc2UgXCJtXCI6XG4gICAgICAgICAgICAgICAgWCArPSBjb21tYW5kWzFdO1xuICAgICAgICAgICAgICAgIFkgKz0gY29tbWFuZFsyXTtcbiAgICAgICAgICAgICAgICBwcmV2UG9pbnQgPSB1bmRlZmluZWQ7XG4gICAgICAgICAgICAgICAgYmVnaW5pbmdQYXRoID0gbmJQb2ludHM7XG4gICAgICAgICAgICAgICAgYnJlYWs7XG4gICAgICAgICAgICBjYXNlIFwiTVwiOlxuICAgICAgICAgICAgICAgIFggPSBjb21tYW5kWzFdO1xuICAgICAgICAgICAgICAgIFkgPSBjb21tYW5kWzJdO1xuICAgICAgICAgICAgICAgIHByZXZQb2ludCA9IHVuZGVmaW5lZDtcbiAgICAgICAgICAgICAgICBiZWdpbmluZ1BhdGggPSBuYlBvaW50cztcbiAgICAgICAgICAgICAgICBicmVhazsgXG4gICAgICAgICAgICBjYXNlIFwiY1wiOlxuICAgICAgICAgICAgICAgIFggPSBjb21tYW5kWzFdO1xuICAgICAgICAgICAgICAgIFkgPSBjb21tYW5kWzJdO1xuICAgICAgICAgICAgICAgIHBvaW50cy5wdXNoKHtpZDpuYlBvaW50cywgeDpYLCB5Oll9KTtcblxuICAgICAgICAgICAgICAgIGlmIChwcmV2UG9pbnQgIT0gdW5kZWZpbmVkKSB7XG4gICAgICAgICAgICAgICAgICAgIGVkZ2VzW3ByZXZQb2ludF0gPSBuYlBvaW50cztcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgcHJldlBvaW50ID0gbmJQb2ludHM7XG4gICAgICAgICAgICAgICAgbmJQb2ludHMrKztcbiAgICAgICAgICAgICAgICBicmVhazsgXG4gICAgICAgICAgICBjYXNlIFwibFwiOlxuICAgICAgICAgICAgICAgIFggPSBjb21tYW5kWzFdO1xuICAgICAgICAgICAgICAgIFkgPSBjb21tYW5kWzJdO1xuICAgICAgICAgICAgICAgIHBvaW50cy5wdXNoKHtpZDpuYlBvaW50cywgeDpYLCB5Oll9KTtcblxuICAgICAgICAgICAgICAgIGlmIChwcmV2UG9pbnQgIT0gdW5kZWZpbmVkKSB7XG4gICAgICAgICAgICAgICAgICAgIGVkZ2VzW3ByZXZQb2ludF0gPSBuYlBvaW50cztcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgcHJldlBvaW50ID0gbmJQb2ludHM7XG4gICAgICAgICAgICAgICAgbmJQb2ludHMrKztcbiAgICAgICAgICAgICAgICBicmVhaztcbiAgICAgICAgICAgIGNhc2UgXCJ6XCI6XG4gICAgICAgICAgICAgICAgZWRnZXNbcHJldlBvaW50XSA9IG5iUG9pbnRzO1xuICAgICAgICAgICAgICAgIGJlZ2luaW5nUGF0aCA9IHVuZGVmaW5lZDtcbiAgICAgICAgICAgICAgICBwcmV2UG9pbnQgPSB1bmRlZmluZWQ7XG4gICAgICAgICAgICAgICAgYnJlYWs7ICAgIFxuICAgICAgICB9XG4gICAgfVxuICAgIHJldHVybiB7cG9pbnRzIDogcG9pbnRzLCBlZGdlcyA6IGVkZ2VzfTtcbn1cblxuLy8gaW5pdGlhbGl6ZSBwb2ludHNcblxubW9kdWxlLmV4cG9ydHMgPSBmdW5jdGlvbihuYlN0YXJ0UG9pbnRzLCBuYlJhbmRvbVBvaW50cyl7XG4gICAgdmFyIHBvaW50cyA9IFtdO1xuICAgIHZhciBmb3JjZWRFZGdlcztcbiAgICB2YXIgY2l0eVNldDtcblxuICAgIGlmICh0ZXh0TWVzaCl7XG5cbiAgICAgICAgdmFyIG15VGV4dCA9IHN2Z1RvUG9pbnRzKHN2Z1N0cmluZyk7XG4gICAgICAgIHBvaW50cyA9IG15VGV4dC5wb2ludHM7XG4gICAgICAgIGZvcmNlZEVkZ2VzID0gbXlUZXh0LmVkZ2VzO1xuICAgICAgICBjaXR5U2V0ID0gcmFuZ2UoMCwgcG9pbnRzLmxlbmd0aCk7XG5cbiAgICAgICAgdmFyIHNjYWxlWCA9IDAuNTtcbiAgICAgICAgdmFyIHNjYWxlWSA9IDAuNTtcbiAgICAgICAgdmFyIGRlbHRhWCA9IDAuMjU7XG4gICAgICAgIHZhciBkZWx0YVkgPSAwLjI7XG5cbiAgICAgICAgLy8gc2NhbGUgcG9pbnRzIHRvIFswLDFdICsgZGVsdGFcbiAgICAgICAgdmFyIG1heFggPSBNYXRoLm1heC5hcHBseShNYXRoLCBwb2ludHMubWFwKGZ1bmN0aW9uKHApe3JldHVybiBwLnh9KSk7XG4gICAgICAgIHZhciBtaW5YID0gTWF0aC5taW4uYXBwbHkoTWF0aCwgcG9pbnRzLm1hcChmdW5jdGlvbihwKXtyZXR1cm4gcC54fSkpO1xuICAgICAgICB2YXIgbWF4WSA9IE1hdGgubWF4LmFwcGx5KE1hdGgsIHBvaW50cy5tYXAoZnVuY3Rpb24ocCl7cmV0dXJuIHAueX0pKTtcbiAgICAgICAgdmFyIG1pblkgPSBNYXRoLm1pbi5hcHBseShNYXRoLCBwb2ludHMubWFwKGZ1bmN0aW9uKHApe3JldHVybiBwLnl9KSk7XG4gICAgICAgIHBvaW50cyA9IHBvaW50cy5tYXAoZnVuY3Rpb24ocCl7XG4gICAgICAgICAgICB2YXIgeCA9IHNjYWxlWCAqIChwLngtbWluWCkvKG1heFgtbWluWCkgKyBkZWx0YVg7XG4gICAgICAgICAgICB2YXIgeSA9IHNjYWxlWSAqIChwLnktbWluWSkvKG1heFktbWluWSkgKyBkZWx0YVk7XG4gICAgICAgICAgICB2YXIgbmV3UG9pbnQgPSBuZXcgUG9pbnQoeCwgeSk7XG4gICAgICAgICAgICBuZXdQb2ludC5pZCA9IHAuaWQ7XG5cbiAgICAgICAgICAgIHJldHVybiBuZXdQb2ludDtcbiAgICAgICAgfSk7XG5cbiAgICAgICAgLy8gb25seSBhZGQgcmFuZG9tIHBvaW50c1xuICAgICAgICB2YXIgbmJQb2ludHMgPSBwb2ludHMubGVuZ3RoO1xuICAgICAgICBmb3IodmFyIGk9MDsgaTxuYlJhbmRvbVBvaW50czsgKytpKSB7XG5cbiAgICAgICAgICAgIHZhciB4ID0gcmFuZG9tKCk7XG4gICAgICAgICAgICB2YXIgeSA9IHJhbmRvbSgpO1xuXG4gICAgICAgICAgICB2YXIgbmV3UG9pbnQgPSBuZXcgUG9pbnQoeCwgeSk7XG4gICAgICAgICAgICBuZXdQb2ludC5pZCA9IG5iUG9pbnRzO1xuXG4gICAgICAgICAgICBwb2ludHMucHVzaChuZXdQb2ludCk7XG5cbiAgICAgICAgICAgIG5iUG9pbnRzKys7XG4gICAgICAgIH1cblxuICAgIH0gZWxzZSB7XG4gICAgICAgIC8vYWRkIHJhbmRvbSBwb2ludHNcblxuICAgICAgICB2YXIgbmJQb2ludHMgPSAwO1xuICAgICAgICBmb3IodmFyIGk9MDsgaTxuYlJhbmRvbVBvaW50czsgKytpKSB7XG5cbiAgICAgICAgICAgIHZhciB4ID0gcmFuZG9tKCk7XG4gICAgICAgICAgICB2YXIgeSA9IHJhbmRvbSgpO1xuXG4gICAgICAgICAgICB2YXIgbmV3UG9pbnQgPSBuZXcgUG9pbnQoeCwgeSk7XG4gICAgICAgICAgICBuZXdQb2ludC5pZCA9IG5iUG9pbnRzO1xuXG4gICAgICAgICAgICBwb2ludHMucHVzaChuZXdQb2ludCk7XG4gICAgICAgICAgICBcbiAgICAgICAgICAgIG5iUG9pbnRzKys7XG4gICAgICAgIH1cblxuICAgICAgICBjaXR5U2V0ID0gcmFuZ2UoMCwgbmJDaXR5KTtcbiAgICAgICAgY29uc29sZS5sb2coY2l0eVNldCk7XG4gICAgfVxuXG5cbiAgICAvLyBpbml0aWFsaXplIHN0YXJ0IHBvaW50c1xuICAgIHZhciBwb3NzaWJsZVN0YXJ0UG9pbnRzSWQgPSBbXTtcblxuICAgIGZvciAodmFyIGkgPSAwOyBpIDwgbmJTdGFydFBvaW50czsgaSsrKXtcbiAgICAgICAgcG9zc2libGVTdGFydFBvaW50c0lkLnB1c2goTWF0aC5mbG9vcihuYlBvaW50cyAqIHJhbmRvbSgpKSk7XG4gICAgfVxuICAgIFxuXG4gICAgcmV0dXJuIHtcbiAgICAgICAgdGV4dE1lc2g6IHRleHRNZXNoLFxuICAgICAgICBwb2ludHM6IHBvaW50cyxcbiAgICAgICAgY2l0eVNldDogY2l0eVNldCxcbiAgICAgICAgcG9zc2libGVTdGFydFBvaW50c0lkOiBwb3NzaWJsZVN0YXJ0UG9pbnRzSWQsXG4gICAgICAgIG5iUmFuZG9tUG9pbnRzOiBuYlJhbmRvbVBvaW50cyxcbiAgICAgICAgZm9yY2VkRWRnZXM6IGZvcmNlZEVkZ2VzXG4gICAgfTtcbn1cbiIsIid1c2Ugc3RyaWN0J1xuXG5tb2R1bGUuZXhwb3J0cyA9IGZ1bmN0aW9uIChjb250YWluZXIpe1xuXG5cdHZhciBtb3VzZSA9IHtcblx0ICAgIHg6IDAsXG5cdCAgICB5OiAwXG5cdH07XG5cblx0Y29udGFpbmVyLmFkZEV2ZW50TGlzdGVuZXIoICdtb3VzZW1vdmUnLCBmdW5jdGlvbihlKXtcblx0ICAgIHZhciByZWN0ID0gY29udGFpbmVyLmdldEJvdW5kaW5nQ2xpZW50UmVjdCgpO1xuXG5cdCAgICBtb3VzZS54ID0gKGUuY2xpZW50WCAtIHJlY3QubGVmdCApIC8gcmVjdC53aWR0aDtcblx0ICAgIG1vdXNlLnkgPSAoZS5jbGllbnRZIC0gcmVjdC50b3AgKS8gcmVjdC5oZWlnaHQ7XG5cdH0pO1xuXG5cdHJldHVybiBtb3VzZTtcblxufTtcbiIsIid1c2Ugc3RyaWN0J1xuXG5mdW5jdGlvbiBQb2ludCh4LCB5KSB7XG4gICAgdGhpcy5pZCA9IHVuZGVmaW5lZDsgICAgICAgICAgICAgICAgXG4gICAgdGhpcy54ID0geDtcbiAgICB0aGlzLnkgPSB5O1xuICAgIHRoaXMubmV4dHMgPSBbXTtcbn1cblxubW9kdWxlLmV4cG9ydHMgPSBQb2ludDsiLCIndXNlIHN0cmljdCdcblxudmFyIGFudEZ1bmN0aW9uID0gcmVxdWlyZSgnLi9hbnQuanMnKTtcbnZhciBhbnRzR3JvdXAgPSByZXF1aXJlKCcuL2FudHNHcm91cCcpO1xuXG52YXIgcmFuZG9tID0gTWF0aC5yYW5kb207XG5cbnZhciBSQU5ET01NVlQgPSAwLjAwMztcbnZhciBBTlRTSVpFID0gMC4wMDI7XG5cbm1vZHVsZS5leHBvcnRzID0gZnVuY3Rpb24oY29udGFpbmVyLCBwb2ludHNNYXAsIG9wdGlvbnMpe1xuXG5cdGlmKCFjb250YWluZXIpXG5cdFx0dGhyb3cgbmV3IFR5cGVFcnJvcignTWlzc2luZyBjb250YWluZXInKTtcblxuXHQvLyBBbnRzIHZhcmlhYmxlc1xuXHR2YXIgZWRnZXMgPSBwb2ludHNNYXAuZWRnZXM7XG5cdHZhciBvYmpQb3B1bGF0aW9uSW5pdGlhbCA9IG9wdGlvbnMubmJBbnRzO1xuXHR2YXIgb2JqUG9wdWxhdGlvbiA9IG9ialBvcHVsYXRpb25Jbml0aWFsO1xuXHR2YXIgcG9pbnRzSW5mb3MgPSBwb2ludHNNYXAucG9pbnRzSW5mb3M7XG5cdHZhciBwb3B1bGF0aW9uID0gW107XG5cdHZhciBuYkFudHNQZXJTdGVwID0gMTAwO1xuXHRcblx0dmFyIEFudCA9IGFudEZ1bmN0aW9uKGNvbnRhaW5lciwgcG9pbnRzSW5mb3MsIG9wdGlvbnMpO1xuXHRhbnRzR3JvdXAgPSBhbnRzR3JvdXAoQW50KTtcblxuXHQvLyBBbmltYXRpb24gdmFyaWFibGVzXG5cdHZhciBhbmltSUQ7XG5cdHZhciBkZWx0YVRpbWU7XG5cdHZhciBGUFNDb3VudDtcblx0dmFyIGxhc3RVcGRhdGUgPSBwZXJmb3JtYW5jZS5ub3coKTtcblx0Ly8gdmFyIEZQU01vbml0b3IgPSBkb2N1bWVudC5xdWVyeVNlbGVjdG9yKCcjRlBTJyk7XG5cdC8vIHZhciBkVE1vbml0b3IgPSBkb2N1bWVudC5xdWVyeVNlbGVjdG9yKCcjZFQnKTtcblx0dmFyIHJlZnJlc2hUaW1lID0gMDtcblx0dmFyIG1heERlbHRhVGltZSA9IDQwO1xuXHR2YXIgRlBTT3ZlckxpbWl0Q291bnQgPSAwO1xuXHR2YXIgRlBTVW5kZXJMaW1pdENvdW50ID0gMDtcblxuXG5cdC8vIENhbnZhc1xuXHR2YXIgY2FudmFzTGlzdCA9IGRvY3VtZW50LmdldEVsZW1lbnRzQnlUYWdOYW1lKFwiY2FudmFzXCIpO1xuXHRcblx0aWYgKGNhbnZhc0xpc3QubGVuZ3RoID09PSAwKXtcblx0XHR2YXIgY2FudmFzID0gZG9jdW1lbnQuY3JlYXRlRWxlbWVudChcImNhbnZhc1wiKTtcblx0XHR2YXIgcmVjdCA9IGNvbnRhaW5lci5nZXRCb3VuZGluZ0NsaWVudFJlY3QoKTtcblx0XHRjYW52YXMud2lkdGggPSByZWN0LndpZHRoO1xuXHRcdGNhbnZhcy5oZWlnaHQgPSByZWN0LmhlaWdodDtcblx0XHRjYW52YXMuc3R5bGUuYmFja2dyb3VuZENvbG9yID0gXCJyZ2JhKDI1MCwgMjUwLCAyNTAsIDApXCI7IFxuXHRcdGNvbnRhaW5lci5hcHBlbmRDaGlsZChjYW52YXMpO1xuXHR9XG5cdGVsc2V7XG5cdFx0dmFyIGNhbnZhcyA9IGNhbnZhc0xpc3RbMF07XG5cdFx0Y29uc29sZS5sb2coJ0NBTlZBUycpO1xuXHR9XG5cdHZhciBjb250ZXh0ID0gY2FudmFzLmdldENvbnRleHQoXCIyZFwiKTtcblx0Y29udGV4dC5jbGVhclJlY3QgKCAwICwgMCAsIGNhbnZhcy53aWR0aCwgY2FudmFzLmhlaWdodCApO1xuXHRcblxuXHRmdW5jdGlvbiBjaGVja0FudE51bWJlcihhbnROdW1iZXIpe1xuXHRcdGlmIChhbnROdW1iZXIgPCBvYmpQb3B1bGF0aW9uIC0gNTApe1xuXHRcdFx0Ly8gRlBTTW9uaXRvci5zdHlsZS5jb2xvciA9IFwiZ3JlZW5cIjtcblx0XHRcdHBvcHVsYXRpb24gPSBhbnRzR3JvdXAuY3JlYXRlKHBvcHVsYXRpb24pO1xuXHRcdH1cdFxuXHRcdGVsc2UgaWYgKGFudE51bWJlciA+IG9ialBvcHVsYXRpb24pe1xuXHRcdFx0cG9wdWxhdGlvbiA9IGFudHNHcm91cC5yZW1vdmUocG9wdWxhdGlvbiwgYW50TnVtYmVyIC0gb2JqUG9wdWxhdGlvbik7XG5cdFx0XHQvLyBGUFNNb25pdG9yLnN0eWxlLmNvbG9yID0gXCJyZWRcIjtcblx0XHR9XG5cdFx0Ly8gZWxzZVxuXHRcdC8vIFx0RlBTTW9uaXRvci5zdHlsZS5jb2xvciA9IFwid2hpdGVcIjtcblx0fVxuXG5cdGZ1bmN0aW9uIGRpc3BsYXlGUFMoZFQpe1xuXHRcdEZQU0NvdW50ID0gKDEwMDAvZFQpLnRvRml4ZWQoMik7XG5cdFx0dmFyIHQgPSBkVC50b0ZpeGVkKDIpO1xuXHRcdC8vIEZQU01vbml0b3IudGV4dENvbnRlbnQgPSAnRlBTIDogJyArIEZQU0NvdW50OyAgXG5cdFx0Ly8gZFRNb25pdG9yLnRleHRDb250ZW50ID0gJ25iQW50cyA6ICcgKyBwb3B1bGF0aW9uLmxlbmd0aDtcblx0XHQvLyBkVE1vbml0b3IuaW5uZXJUZXh0ID0gJ2RUIDogJyArIHQgKyAnbXMnO1xuXHR9XG5cblx0ZnVuY3Rpb24gdGljaygpIHtcblx0XHR2YXIgbm93ID0gcGVyZm9ybWFuY2Uubm93KCk7XG5cdFx0ZGVsdGFUaW1lID0gbm93IC0gbGFzdFVwZGF0ZTtcblx0XHRsYXN0VXBkYXRlID0gbm93O1xuXHRcdHJlZnJlc2hUaW1lICs9IGRlbHRhVGltZS8xMDAwOyAvLyBpbiBzZWNvbmRzXG5cblx0XHQvLyBjb25zb2xlLmxvZygnbmJBbnRzJywgcG9wdWxhdGlvbi5sZW5ndGgpO1xuXG5cdFx0Y2hlY2tBbnROdW1iZXIocG9wdWxhdGlvbi5sZW5ndGgpO1xuXG5cdFx0Ly8gZGlzcGxheSBGUFMgaW5mbyBldmVyeSAwLjMgc1xuXHRcdGlmIChyZWZyZXNoVGltZSA+IDAuMyl7XG5cdFx0XHRkaXNwbGF5RlBTKGRlbHRhVGltZSk7XG5cdFx0XHRyZWZyZXNoVGltZSA9IDA7IFxuXHRcdH1cblxuXHRcdC8vIHJlbW92ZSBhbnRzIHdoZW4gZnJhbWUgcmF0ZSBpcyB0b28gbG93XG5cdFx0aWYgKEZQU092ZXJMaW1pdENvdW50ID09PSAxMCkge1xuXHRcdFx0b2JqUG9wdWxhdGlvbiA9IG9ialBvcHVsYXRpb24gKiBtYXhEZWx0YVRpbWUgLyBkZWx0YVRpbWU7XG5cdFx0XHRGUFNPdmVyTGltaXRDb3VudCA9IDA7XG5cdFx0fVxuXG5cdFx0d2hpbGUgKEZQU1VuZGVyTGltaXRDb3VudCA+IDUwICYmIG9ialBvcHVsYXRpb24gPCBvYmpQb3B1bGF0aW9uSW5pdGlhbCkge1xuXHRcdFx0b2JqUG9wdWxhdGlvbiArPSAxMDtcblx0XHR9XG5cblx0XHQvLyBjaGVjayBkdXJhdGlvbiBvZiBvdmVyL3VuZGVyIGZyYW1lcmF0ZSBsaW1pdCBwZXJpb2RzXG5cdFx0aWYgKGRlbHRhVGltZSA+IG1heERlbHRhVGltZSl7XG5cdFx0XHRGUFNPdmVyTGltaXRDb3VudCsrO1xuXHRcdFx0RlBTVW5kZXJMaW1pdENvdW50ID0gMDtcblx0XHR9XG5cdFx0ZWxzZSB7XG5cdFx0XHRGUFNPdmVyTGltaXRDb3VudCA9IDA7XG5cdFx0XHRGUFNVbmRlckxpbWl0Q291bnQrKztcblx0XHR9XG5cblx0XHQvLyBkcmF3IGluIGNhbnZhc1xuXHRcdHZhciB3ID0gY2FudmFzLndpZHRoO1xuXHRcdHZhciBoID0gY2FudmFzLmhlaWdodDtcblx0XHR2YXIgbW91c2UgPSBbbGFzdE1vdXNlTW92ZUV2ZW50LmNsaWVudFgvdywgbGFzdE1vdXNlTW92ZUV2ZW50LmNsaWVudFkvaF07XG5cdFx0Y29udGV4dC5zZXRUcmFuc2Zvcm0odywgMCwgMCwgaCwgMCwgMCk7XG5cdFx0Y29udGV4dC5maWxsU3R5bGUgPSBcInJnYmEoMjUwLCAyNTAsIDI1MCwgMC40KVwiO1xuXHRcdGNvbnRleHQuZmlsbFJlY3QoMCwwLHcsaCk7XG5cblx0XHQvLyAvLyBlZGdlc1xuXHRcdC8vIGNvbnRleHQuc3Ryb2tlU3R5bGUgPSBcIiMwMDBcIjtcblx0XHQvLyBmb3IodmFyIGk9MDsgaSA8IGVkZ2VzLmxlbmd0aDsgKytpKSB7XG5cdFx0Ly8gICAgIGNvbnRleHQubGluZVdpZHRoID0gMC4wMDAxO1xuXHRcdC8vICAgICB2YXIgZWRnZSA9IGVkZ2VzW2ldO1xuXHRcdC8vICAgICBpZiAoZWRnZS5waGVyb21vbiAhPSAwKXtcblx0XHQvLyAgICAgICAgIGNvbnRleHQubGluZVdpZHRoID0gTWF0aC5taW4oMC4wMDAwMSAqIGVkZ2UucGhlcm9tb24sIDAuMDEpO1xuXHRcdC8vICAgICB9IGVsc2Uge1xuXHRcdC8vICAgICAgICAgY29udGV4dC5saW5lV2lkdGggPSAwLjAwMDAxO1xuXHRcdC8vICAgICB9XG5cdFx0Ly8gICAgIGNvbnRleHQuYmVnaW5QYXRoKCk7XG5cdFx0Ly8gICAgIGNvbnRleHQubW92ZVRvKHBvaW50c0luZm9zLnBvaW50c1tlZGdlLnB0MS5pZF0ueCwgcG9pbnRzSW5mb3MucG9pbnRzW2VkZ2UucHQxLmlkXS55KTtcblx0XHQvLyAgICAgY29udGV4dC5saW5lVG8ocG9pbnRzSW5mb3MucG9pbnRzW2VkZ2UucHQyLmlkXS54LCBwb2ludHNJbmZvcy5wb2ludHNbZWRnZS5wdDIuaWRdLnkpO1xuXHRcdC8vICAgICBjb250ZXh0LnN0cm9rZSgpO1xuXHRcdC8vIH1cblxuXHRcdC8vIC8vIHZlcnRpY2VzXG5cdFx0Ly8gZm9yKHZhciBpPTA7IGk8cG9pbnRzSW5mb3MucG9pbnRzLmxlbmd0aDsgKytpKSB7XG5cdFx0Ly8gICAgIGNvbnRleHQuYmVnaW5QYXRoKClcblx0XHQvLyAgICAgdmFyIHBvaW50ID0gcG9pbnRzSW5mb3MucG9pbnRzW2ldO1xuXHRcdC8vICAgICBpZiAocG9pbnRzSW5mb3MuY2l0eVNldC5pbmRleE9mKHBvaW50LmlkKSAhPSAtMSl7XG5cdFx0Ly8gICAgICAgICBjb250ZXh0LmZpbGxTdHlsZSA9IFwiIzAxMDFERlwiO1xuXHRcdC8vICAgICAgICAgY29udGV4dC5hcmMocG9pbnQueCwgcG9pbnQueSwgMC4wMDYsIDAsIDIqTWF0aC5QSSk7XG5cdFx0Ly8gICAgIH1cblx0XHQvLyAgICAgZWxzZSB7XG5cdFx0Ly8gICAgICAgICBjb250ZXh0LmZpbGxTdHlsZSA9IFwiIzAwMFwiO1xuXHRcdC8vICAgICAgICAgY29udGV4dC5hcmMocG9pbnRzSW5mb3MucG9pbnRzW2ldLngsIHBvaW50c0luZm9zLnBvaW50c1tpXS55LCAwLjAwMywgMCwgMipNYXRoLlBJKTtcblx0XHQvLyAgICAgfVxuXHRcdC8vICAgICBjb250ZXh0LmNsb3NlUGF0aCgpO1xuXHRcdC8vICAgICBjb250ZXh0LmZpbGwoKTtcblx0XHQvLyB9XG5cblx0XHQvLyBtb3ZlIGFudHNcblx0XHRwb3B1bGF0aW9uLmZvckVhY2goZnVuY3Rpb24oYW50KXtcblx0XHRcdGFudC50cmFuc2l0KCk7XG5cdFx0fSk7XG5cblx0XHQvLyBwaGVyb21vbiBldmFwb3JhdGlvblxuXHRcdGVkZ2VzLmZvckVhY2goZnVuY3Rpb24oZWRnZSl7XG5cdFx0XHRpZihlZGdlLnBoZXJvbW9uID4gMCl7XG5cdFx0XHRcdGVkZ2UucGhlcm9tb24gLT0gMC4wMDAxO1xuXHRcdFx0fVxuXHRcdH0pO1xuXG5cdFx0Ly8gYW50c1xuXHRcdHBvcHVsYXRpb24uZm9yRWFjaChmdW5jdGlvbihhbnQpe1xuXHRcdFx0Y29udGV4dC5iZWdpblBhdGgoKVxuXHRcdFx0dmFyIHggPSBhbnQueCArIFJBTkRPTU1WVCpyYW5kb20oKTtcblx0XHRcdHZhciB5ID0gYW50LnkgKyBSQU5ET01NVlQqcmFuZG9tKCk7XG5cblx0XHRcdGNvbnRleHQuZmlsbFN0eWxlID0gXCJibGFja1wiXG5cdFx0XHRjb250ZXh0LmZpbGxSZWN0KHgsIHksIEFOVFNJWkUsIEFOVFNJWkUpO1xuXHRcdFx0Y29udGV4dC5jbG9zZVBhdGgoKTtcblx0XHRcdGNvbnRleHQuZmlsbCgpO1xuXHRcdH0pXG5cdH07XG5cdFxuXHR2YXIgbGFzdE1vdXNlTW92ZUV2ZW50ID0ge1xuXHRcdGNsaWVudFg6IDAsXG5cdFx0Y2xpZW50WTogMFxuXHR9O1xuXHRcblx0Y29udGFpbmVyLmFkZEV2ZW50TGlzdGVuZXIoJ21vdXNlbW92ZScsIGZ1bmN0aW9uKGUpe1xuXHRcdGxhc3RNb3VzZU1vdmVFdmVudCA9IGU7XG5cdH0pO1xuXHRcblx0dmFyIHBhdXNlZCA9IGZhbHNlO1xuXHRcblx0ZnVuY3Rpb24gdG9nZ2xlUGxheVBhdXNlKCl7XG5cdFx0cGF1c2VkID0gIXBhdXNlZDtcblx0XHRpZighcGF1c2VkKVxuXHRcdFx0YW5pbWF0ZSgpO1xuXHR9XG5cblx0ZnVuY3Rpb24gcmVzZXQoKXtcblx0XHRwb3B1bGF0aW9uID0gW107XG5cdFx0ZWRnZXMgPSBbXTtcblx0XHRwb2ludHNJbmZvcyA9IFtdO1xuXG5cdFx0Y2FuY2VsQW5pbWF0aW9uRnJhbWUoYW5pbUlEKTtcblx0fVxuXHRcblx0Ly8gY29udGFpbmVyLmFkZEV2ZW50TGlzdGVuZXIoJ2NsaWNrJywgdG9nZ2xlUGxheVBhdXNlKTtcblxuXHRmdW5jdGlvbiBhbmltYXRlKCl7XG5cdFx0dGljaygpO1xuXHRcdFxuXHRcdGlmKCFwYXVzZWQpXG5cdFx0XHRhbmltSUQgPSByZXF1ZXN0QW5pbWF0aW9uRnJhbWUoYW5pbWF0ZSk7XG5cdH1cblx0YW5pbWF0ZSgpO1xuXG5cdGZ1bmN0aW9uIG1vZGlmeUFudHMob3B0cyl7XG5cdFx0b2JqUG9wdWxhdGlvbiA9IG9wdHMubmJBbnRzO1xuXG5cdFx0cG9wdWxhdGlvbi5mb3JFYWNoKGZ1bmN0aW9uKGFudCl7XG5cdFx0XHRhbnQudmVsb2NpdHkgPSBvcHRzLnZlbG9jaXR5O1xuXHRcdFx0YW50LndlaWdodCA9IG9wdHMud2VpZ2h0O1xuXHRcdFx0YW50LnJlcFNpemUgPSBvcHRzLnJlcFNpemU7XG5cdFx0XHRhbnQucmVwU3BlZWQgPSBvcHRzLnJlcFNwZWVkO1xuXHRcdH0pO1xuXHR9XG5cdFxuXHRyZXR1cm4ge1xuXHRcdHRvZ2dsZVBsYXlQYXVzZTogdG9nZ2xlUGxheVBhdXNlLFxuXHRcdHJlc2V0OiByZXNldCxcblx0XHQvLyBzaG91bGQgYmUgYSBnZXR0ZXIvc2V0dGVyLCBidXQgSUU4XG5cdFx0Z2V0QW50Q291bnQ6IGZ1bmN0aW9uKCl7XG5cdFx0XHRyZXR1cm4gcG9wdWxhdGlvbi5sZW5ndGg7XG5cdFx0fSxcblx0XHRtb2RpZnlBbnRzOiBtb2RpZnlBbnRzXG5cdH1cbn1cbiIsIid1c2Ugc3RyaWN0JztcblxudmFyIHBvaW50MSA9IFwibSAxOC4yNSwxOS41IGMgMTguMjUsMTUuNzY3NzY4OCAxNS4wMDAwOTUxLDEyLjc1IDExLDEyLjc1IGMgNi45OTk5MDQ4OCwxMi43NSAzLjc1LDE1Ljc2Nzc2ODggMy43NSwxOS41IGMgMy43NSwyMy4yMzIyMzEyIDYuOTk5OTA0ODgsMjYuMjUgMTEsMjYuMjUgYyAxNS4wMDAwOTUxLDI2LjI1IDE4LjI1LDIzLjIzMjIzMTIgMTguMjUsMTkuNSB6IG0gNC4yNSwxOS41IGMgNC4yNSwxNi4wNTI1Mjk0IDcuMjY4MTA4NjMsMTMuMjUgMTEsMTMuMjUgYyAxNC43MzE4OTE0LDEzLjI1IDE3Ljc1LDE2LjA1MjUyOTQgMTcuNzUsMTkuNSBjIDE3Ljc1LDIyLjk0NzQ3MDYgMTQuNzMxODkxNCwyNS43NSAxMSwyNS43NSBjIDcuMjY4MTA4NjMsMjUuNzUgNC4yNSwyMi45NDc0NzA2IDQuMjUsMTkuNSB6XCI7XG52YXIgcG9pbnQyID0gXCJtIDg5LjI1LDguNSBjIDg5LjI1LDQuMjE2MDU3ODcgODUuNTUyODcxNiwwLjc1IDgxLDAuNzUgYyA3Ni40NDcxMjg0LDAuNzUgNzIuNzUsNC4yMTYwNTc4NyA3Mi43NSw4LjUgYyA3Mi43NSwxMi43ODM5NDIxIDc2LjQ0NzEyODQsMTYuMjUgODEsMTYuMjUgYyA4NS41NTI4NzE2LDE2LjI1IDg5LjI1LDEyLjc4Mzk0MjEgODkuMjUsOC41IHogbSA3My4yNSw4LjUgYyA3My4yNSw0LjQ5OTY3MDg4IDc2LjcxNjMxNTYsMS4yNSA4MSwxLjI1IGMgODUuMjgzNjg0NCwxLjI1IDg4Ljc1LDQuNDk5NjcwODggODguNzUsOC41IGMgODguNzUsMTIuNTAwMzI5MSA4NS4yODM2ODQ0LDE1Ljc1IDgxLDE1Ljc1IGMgNzYuNzE2MzE1NiwxNS43NSA3My4yNSwxMi41MDAzMjkxIDczLjI1LDguNSB6XCI7XG52YXIgcG9pbnQzID0gXCJtIDE2MC4yNSwxMSBjIDE2MC4yNSw1LjMzMzQ4Mzc1IDE1NS4yMDgxNjgsMC43NSAxNDksMC43NSBjIDE0Mi43OTE4MzIsMC43NSAxMzcuNzUsNS4zMzM0ODM3NSAxMzcuNzUsMTEgYyAxMzcuNzUsMTYuNjY2NTE2MiAxNDIuNzkxODMyLDIxLjI1IDE0OSwyMS4yNSBjIDE1NS4yMDgxNjgsMjEuMjUgMTYwLjI1LDE2LjY2NjUxNjIgMTYwLjI1LDExIHogbSAxMzguMjUsMTEgYyAxMzguMjUsNS42MjA4MjEyNSAxNDMuMDU3OTAzLDEuMjUgMTQ5LDEuMjUgYyAxNTQuOTQyMDk3LDEuMjUgMTU5Ljc1LDUuNjIwODIxMjUgMTU5Ljc1LDExIGMgMTU5Ljc1LDE2LjM3OTE3ODcgMTU0Ljk0MjA5NywyMC43NSAxNDksMjAuNzUgYyAxNDMuMDU3OTAzLDIwLjc1IDEzOC4yNSwxNi4zNzkxNzg3IDEzOC4yNSwxMSB6XCI7XG52YXIgcG9pbnQ0ID0gXCJtIDE2MC4yNSw3Ni41IGMgMTYwLjI1LDcyLjc2Nzc2ODggMTU3LjAwMDA5NSw2OS43NSAxNTMsNjkuNzUgYyAxNDguOTk5OTA1LDY5Ljc1IDE0NS43NSw3Mi43Njc3Njg4IDE0NS43NSw3Ni41IGMgMTQ1Ljc1LDgwLjIzMjIzMTIgMTQ4Ljk5OTkwNSw4My4yNSAxNTMsODMuMjUgYyAxNTcuMDAwMDk1LDgzLjI1IDE2MC4yNSw4MC4yMzIyMzEyIDE2MC4yNSw3Ni41IHogbSAxNDYuMjUsNzYuNSBjIDE0Ni4yNSw3My4wNTI1Mjk0IDE0OS4yNjgxMDksNzAuMjUgMTUzLDcwLjI1IGMgMTU2LjczMTg5MSw3MC4yNSAxNTkuNzUsNzMuMDUyNTI5NCAxNTkuNzUsNzYuNSBjIDE1OS43NSw3OS45NDc0NzA2IDE1Ni43MzE4OTEsODIuNzUgMTUzLDgyLjc1IGMgMTQ5LjI2ODEwOSw4Mi43NSAxNDYuMjUsNzkuOTQ3NDcwNiAxNDYuMjUsNzYuNSB6XCI7XG52YXIgcG9pbnQ1ID0gXCJtIDk1LjI1LDc2IGMgOTUuMjUsNzAuMzM2Mjc5NSA5MC40MzQ0MDY1LDY1Ljc1IDg0LjUsNjUuNzUgYyA3OC41NjU1OTM1LDY1Ljc1IDczLjc1LDcwLjMzNjI3OTUgNzMuNzUsNzYgYyA3My43NSw4MS42NjM3MjA1IDc4LjU2NTU5MzUsODYuMjUgODQuNSw4Ni4yNSBjIDkwLjQzNDQwNjUsODYuMjUgOTUuMjUsODEuNjYzNzIwNSA5NS4yNSw3NiB6IG0gNzQuMjUsNzYgYyA3NC4yNSw3MC42MTgwMjU1IDc4LjgzNjQyNjgsNjYuMjUgODQuNSw2Ni4yNSBjIDkwLjE2MzU3MzIsNjYuMjUgOTQuNzUsNzAuNjE4MDI1NSA5NC43NSw3NiBjIDk0Ljc1LDgxLjM4MTk3NDUgOTAuMTYzNTczMiw4NS43NSA4NC41LDg1Ljc1IGMgNzguODM2NDI2OCw4NS43NSA3NC4yNSw4MS4zODE5NzQ1IDc0LjI1LDc2IHpcIjtcbnZhciBwb2ludDYgPSBcIm0gMjAuMjUsNzUgYyAyMC4yNSw3MC45OTE5MzM4IDE2Ljc3NjQ5OTUsNjcuNzUgMTIuNSw2Ny43NSBjIDguMjIzNTAwNDYsNjcuNzUgNC43NSw3MC45OTE5MzM4IDQuNzUsNzUgYyA0Ljc1LDc5LjAwODA2NjIgOC4yMjM1MDA0Niw4Mi4yNSAxMi41LDgyLjI1IGMgMTYuNzc2NDk5NSw4Mi4yNSAyMC4yNSw3OS4wMDgwNjYyIDIwLjI1LDc1IHogbSA1LjI1LDc1IGMgNS4yNSw3MS4yNzYwNzk3IDguNDkyMjI4MjksNjguMjUgMTIuNSw2OC4yNSBjIDE2LjUwNzc3MTcsNjguMjUgMTkuNzUsNzEuMjc2MDc5NyAxOS43NSw3NSBjIDE5Ljc1LDc4LjcyMzkyMDMgMTYuNTA3NzcxNyw4MS43NSAxMi41LDgxLjc1IGMgOC40OTIyMjgyOSw4MS43NSA1LjI1LDc4LjcyMzkyMDMgNS4yNSw3NSB6XCI7XG52YXIgcG9pbnQ3ID0gXCJtIDIzLjI1LDEzOSBjIDIzLjI1LDEzMi43ODY3OTcgMTguMjEzMjAzNCwxMjcuNzUgMTIsMTI3Ljc1IGMgNS43ODY3OTY1NiwxMjcuNzUgMC43NSwxMzIuNzg2Nzk3IDAuNzUsMTM5IGMgMC43NSwxNDUuMjEzMjAzIDUuNzg2Nzk2NTYsMTUwLjI1IDEyLDE1MC4yNSBjIDE4LjIxMzIwMzQsMTUwLjI1IDIzLjI1LDE0NS4yMTMyMDMgMjMuMjUsMTM5IHogbSAxLjI1LDEzOSBjIDEuMjUsMTMzLjA2MjkzOSA2LjA2MjkzODk0LDEyOC4yNSAxMiwxMjguMjUgYyAxNy45MzcwNjExLDEyOC4yNSAyMi43NSwxMzMuMDYyOTM5IDIyLjc1LDEzOSBjIDIyLjc1LDE0NC45MzcwNjEgMTcuOTM3MDYxMSwxNDkuNzUgMTIsMTQ5Ljc1IGMgNi4wNjI5Mzg5NCwxNDkuNzUgMS4yNSwxNDQuOTM3MDYxIDEuMjUsMTM5IHpcIjtcbnZhciBwb2ludDggPSBcIm0gOTUuMjQyOTIwOSwxMzMuMTQxNjI0IGMgOTUuNTUzNzkyNywxMjcuMjA5ODM2IDkwLjc1NjI2OTgsMTIyLjE0MTczNiA4NC41MzM3NDA0LDEyMS44MTU2MjcgYyA3OC4zMTEyMTEsMTIxLjQ4OTUxOCA3My4wMTAyMDg2LDEyNi4wMjgzNzcgNzIuNjk5MzM2OCwxMzEuOTYwMTY1IGMgNzIuMzg4NDY1LDEzNy44OTE5NTMgNzcuMTg1OTg3OSwxNDIuOTYwMDUzIDgzLjQwODUxNzMsMTQzLjI4NjE2MiBjIDg5LjYzMTA0NjcsMTQzLjYxMjI3MSA5NC45MzIwNDksMTM5LjA3MzQxMiA5NS4yNDI5MjA5LDEzMy4xNDE2MjQgeiBtIDczLjE5ODY1MTYsMTMxLjk4NjMzMyBjIDczLjQ5NDc3MTEsMTI2LjMzNjAzNiA3OC41NTUzODgsMTIyLjAwMzAwMSA4NC41MDc1NzI0LDEyMi4zMTQ5NDIgYyA5MC40NTk3NTY3LDEyMi42MjY4ODMgOTUuMDM5NzI1NiwxMjcuNDY1MTU5IDk0Ljc0MzYwNjEsMTMzLjExNTQ1NiBjIDk0LjQ0NzQ4NjYsMTM4Ljc2NTc1NCA4OS4zODY4Njk2LDE0My4wOTg3ODggODMuNDM0Njg1MywxNDIuNzg2ODQ3IGMgNzcuNDgyNTAwOSwxNDIuNDc0OTA3IDcyLjkwMjUzMiwxMzcuNjM2NjMgNzMuMTk4NjUxNiwxMzEuOTg2MzMzIHpcIjtcbnZhciBwb2ludDkgPSBcIm0gMTY3LjcyODAyMywxMzUuMDYyNDc5IGMgMTY4LjAzODMyNywxMjkuMTQxNTMzIDE2My40ODIxMTcsMTI0LjA4OTg2NCAxNTcuNTUxNjYyLDEyMy43NzkwNjIgYyAxNTEuNjIxMjA4LDEyMy40NjgyNjEgMTQ2LjU2MTkxNCwxMjguMDE2MDAyIDE0Ni4yNTE2MSwxMzMuOTM2OTQ4IGMgMTQ1Ljk0MTMwNywxMzkuODU3ODk0IDE1MC40OTc1MTcsMTQ0LjkwOTU2MiAxNTYuNDI3OTcxLDE0NS4yMjAzNjQgYyAxNjIuMzU4NDI2LDE0NS41MzExNjYgMTY3LjQxNzcyLDE0MC45ODM0MjUgMTY3LjcyODAyMywxMzUuMDYyNDc5IHogbSAxNDYuNzUwOTI1LDEzMy45NjMxMTYgYyAxNDcuMDQ2NzY3LDEyOC4zMTgxMjEgMTUxLjg3MDYxNywxMjMuOTgyMDE4IDE1Ny41MjU0OTQsMTI0LjI3ODM3NyBjIDE2My4xODAzNzIsMTI0LjU3NDczNyAxNjcuNTI0NTUsMTI5LjM5MTMxNyAxNjcuMjI4NzA5LDEzNS4wMzYzMTEgYyAxNjYuOTMyODY3LDE0MC42ODEzMDUgMTYyLjEwOTAxNywxNDUuMDE3NDA5IDE1Ni40NTQxMzksMTQ0LjcyMTA1IGMgMTUwLjc5OTI2MiwxNDQuNDI0NjkgMTQ2LjQ1NTA4MywxMzkuNjA4MTEgMTQ2Ljc1MDkyNSwxMzMuOTYzMTE2IHpcIjtcbnZhciBsZXR0cmVzQU5UID0gXCJtIDI1LjYxODU0MTQsMzEuOTQ0NzI2NiBsIDI1LjYxODU0MTQsMzEuMTAwMzkwNiBsIDI1LjM4ODI2OCwzMC40MDk1NzAzIGwgMjQuOTI3NzIxMSwzMC4wMjU3ODEyIGwgMjQuNTQzOTMyMSwyOS41NjUyMzQ0IGwgMjQuMDgzMzg1MiwyOS4xODE0NDUzIGwgMjMuNjk5NTk2MSwyOC43MjA4OTg0IGwgMjMuMDA4Nzc1OCwyOC40OTA2MjUgbCAyMi4xNDUzNjc4LDI4LjQ5MDYyNSBsIDIxLjE1MTI3MTEsMjguNDkwNjI1IGwgMjAuNDc1NzY4LDI4LjQ5MDYyNSBsIDE5Ljg2MTcwNTUsMjguNzIwODk4NCBsIDE5LjQ3NzkxNjQsMjkuMTgxNDQ1MyBsIDE4Ljc4NzA5NjEsMjkuMzM0OTYwOSBsIDE3LjU0Nzg1NzksMjkuMzM0OTYwOSBsIDE2LjY2OTkwNywyOS4zMzQ5NjA5IGwgMTYuMDAyNDk3MSwyOS4zMzQ5NjA5IGwgMTUuMzM0OTgwMiwyOS4zMzQ5NjA5IGwgMTQuNTY1NDE2NCwyOS4zMzQ5NjA5IGwgMTMuOTUxMzUzOSwyOS41NjUyMzQ0IGwgMTMuNTY3NTY0OSwzMC4wMjU3ODEyIGwgMTIuODc2NzQ0NiwzMC4xNzkyOTY5IGwgMTIuMDMyNDA4NiwzMC4xNzkyOTY5IGwgMTEuNDE4MzQ2MSwzMC40MDk1NzAzIGwgMTEuMDM0NTU3MSwzMC44NzAxMTcyIGwgMTAuMzQzNzM2NywzMS4wMjM2MzI4IGwgOS43Mjk2NzQyNCwzMS4yNTM5MDYyIGwgOS4zNDU4ODUxOCwzMS43MTQ0NTMxIGwgOC44ODUzMzgzLDMyLjA5ODI0MjIgbCA4LjUwMTU0OTI0LDMyLjU1ODc4OTEgbCA4LjA0MTAwMjM3LDMyLjk0MjU3ODEgbCA3Ljg4NzQ4Njc0LDMzLjYzMzM5ODQgbCA3LjY1NzIxMzMsMzQuMjQ3NDYwOSBsIDcuMTk2NjY2NDMsMzQuNjMxMjUgbCA3LjA0MzE1MDgsMzUuMzIyMDcwMyBsIDcuMDQzMTUwOCwzNi4xNjY0MDYyIGwgNy4xOTY2NjY0MywzNi43ODA0Njg3IGwgNy42NTcyMTMzLDM3LjE2NDI1NzggbCA4LjA0MTAwMjM3LDM3LjYyNDgwNDcgbCA4LjUwMTU0OTI0LDM4LjAwODU5MzcgbCA4Ljg4NTMzODMsMzguNDY5MTQwNiBsIDkuMzQ1ODg1MTgsMzguODUyOTI5NyBsIDkuNzI5Njc0MjQsMzkuMzEzNDc2NiBsIDEwLjM0MzczNjcsMzkuNDY2OTkyMiBsIDExLjQ3MjU5NjEsMzkuNDY2OTkyMiBsIDEyLjMxMjI5NjEsMzkuNDY2OTkyMiBsIDEyLjg3NDM4MTgsMzkuNDY2OTkyMiBsIDEzLjY0OTY4NzMsMzkuNDY2OTkyMiBsIDE0LjUzOTg1MjgsMzkuNDY2OTkyMiBsIDE1LjQwOTc1MjQsMzkuNDY2OTkyMiBsIDE2LjEwMDU3MjcsMzkuMzEzNDc2NiBsIDE2LjQ4NDM2MTcsMzguODUyOTI5NyBsIDE3LjA5ODQyNDIsMzguNjIyNjU2MiBsIDE3Ljc4OTI0NDYsMzguNDY5MTQwNiBsIDE4LjE3MzAzMzYsMzguMDA4NTkzNyBsIDE4LjYzMzU4MDUsMzcuNjI0ODA0NyBsIDE5LjAxNzM2OTYsMzcuMTY0MjU3OCBsIDE5LjYzMTQzMjEsMzYuOTMzOTg0NCBsIDIwLjMyMjI1MjQsMzYuNzgwNDY4NyBsIDIwLjcwNjA0MTQsMzYuMzE5OTIxOSBsIDIxLjMyMDEwMzksMzYuMDg5NjQ4NCBsIDIyLjAxMDkyNDIsMzUuOTM2MTMyOCBsIDIyLjM5NDcxMzMsMzUuNDc1NTg1OSBsIDIzLjAwODc3NTgsMzUuMjQ1MzEyNSBsIDIzLjY5OTU5NjEsMzUuMDkxNzk2OSBsIDIzLjkyOTg2OTYsMzQuNDc3NzM0NCBsIDI0LjA4MzM4NTIsMzMuNzg2OTE0MSBsIDI0LjU0MzkzMjEsMzMuNDAzMTI1IGwgMjQuOTI3NzIxMSwzMi45NDI1NzgxIGwgMjUuMzg4MjY4LDMyLjU1ODc4OTEgbCAyNS42MTg1NDE0LDMxLjk0NDcyNjYgeiBtIDM3LjU5Mjc2MDIsNDMuOTE4OTQ1MyBsIDM3LjIwODk3MTEsNDQuMzc5NDkyMiBsIDM2Ljc0ODQyNDIsNDQuNzYzMjgxMyBsIDM2LjM2NDYzNTIsNDUuMjIzODI4MSBsIDM1LjY3MzgxNDksNDUuMzc3MzQzOCBsIDM1LjA1OTc1MjQsNDUuNjA3NjE3MiBsIDM0LjY3NTk2MzMsNDYuMDY4MTY0MSBsIDMzLjk4NTE0Myw0Ni4yMjE2Nzk3IGwgMzMuMzk1MTcyLDQ2LjIyMTY3OTcgbCAzMi40ODc0MjE5LDQ2LjIyMTY3OTcgbCAzMS40NTIxMzUyLDQ2LjIyMTY3OTcgbCAzMC44MzgwNzI3LDQ2LjA2ODE2NDEgbCAzMC40NTQyODM2LDQ1LjYwNzYxNzIgbCAyOS43NjM0NjMzLDQ1LjM3NzM0MzggbCAyOS4xNDk0MDA4LDQ1LjIyMzgyODEgbCAyOC43NjU2MTE3LDQ0Ljc2MzI4MTMgbCAyOC4zMDUwNjQ5LDQ0LjM3OTQ5MjIgbCAyNy45MjEyNzU4LDQzLjkxODk0NTMgbCAyNy4yMzA0NTU1LDQzLjY4ODY3MTkgbCAyNi42MTYzOTMsNDMuNTM1MTU2MyBsIDI2LjIzMjYwMzksNDMuMDc0NjA5NCBsIDI1Ljc3MjA1NzEsNDMuMDc0NjA5NCBsIDI1LjM4ODI2OCw0My41MzUxNTYzIGwgMjQuNjk3NDQ3Nyw0My42ODg2NzE5IGwgMjMuODUzMTExNyw0My42ODg2NzE5IGwgMjMuMjM5MDQ5Miw0My45MTg5NDUzIGwgMjIuODU1MjYwMiw0NC4zNzk0OTIyIGwgMjIuMTY0NDM5OSw0NC41MzMwMDc4IGwgMjEuNTUwMzc3NCw0NC43NjMyODEzIGwgMjEuMTY2NTg4Myw0NS4yMjM4MjgxIGwgMjAuNDc1NzY4LDQ1LjM3NzM0MzggbCAxOS42MzE0MzIxLDQ1LjM3NzM0MzggbCAxOS4wMTczNjk2LDQ1LjYwNzYxNzIgbCAxOC42MzM1ODA1LDQ2LjA2ODE2NDEgbCAxNy45NDI3NjAyLDQ2LjIyMTY3OTcgbCAxNy4xMjc5MDg3LDQ2LjIyMTY3OTcgbCAxNi4yNTQwODgzLDQ2LjIyMTY3OTcgbCAxNS42NDAwMjU4LDQ2LjQ1MTk1MzEgbCAxNS4yNTYyMzY3LDQ2LjkxMjUgbCAxNC43OTU2ODk5LDQ2LjkxMjUgbCAxNC40MTE5MDA4LDQ2LjQ1MTk1MzEgbCAxMy43MjEwODA1LDQ2LjIyMTY3OTcgbCAxMi43OTI2ODM3LDQ2LjIyMTY3OTcgbCAxMi4wMTY5OCw0Ni4yMjE2Nzk3IGwgMTEuMzY0NjksNDYuMjIxNjc5NyBsIDEwLjQ4NDkxLDQ2LjIyMTY3OTcgbCA5LjQ5OTQwMDgsNDYuMjIxNjc5NyBsIDguODg1MzM4Myw0Ni4wNjgxNjQxIGwgOC41MDE1NDkyNCw0NS42MDc2MTcyIGwgNy44MTA3Mjg5Myw0NS4zNzczNDM4IGwgNy4xOTY2NjY0Myw0NS4yMjM4MjgxIGwgNi44MTI4NzczNyw0NC43NjMyODEzIGwgNi4xMjIwNTcwNSw0NC41MzMwMDc4IGwgNS41MDc5OTQ1NSw0NC4zNzk0OTIyIGwgNS4xMjQyMDU0OSw0My45MTg5NDUzIGwgNC42NjM2NTg2Miw0My41MzUxNTYzIGwgNC4yNzk4Njk1NSw0My4wNzQ2MDk0IGwgMy44MTkzMjI2OCw0Mi42OTA4MjAzIGwgMy40MzU1MzM2Miw0Mi4yMzAyNzM0IGwgMi45NzQ5ODY3NCw0MS44NDY0ODQ0IGwgMi41OTExOTc2OCw0MS4zODU5Mzc1IGwgMi4xMzA2NTA4LDQxLjAwMjE0ODQgbCAxLjk3NzEzNTE4LDQwLjM4ODA4NTkgbCAxLjc0Njg2MTc0LDM5LjY5NzI2NTYgbCAxLjI4NjMxNDg3LDM5LjMxMzQ3NjYgbCAxLjEzMjc5OTI0LDM4LjY5OTQxNDEgbCAxLjEzMjc5OTI0LDM3Ljg1NTA3ODEgbCAwLjkwMjUyNTgwNCwzNy4xNjQyNTc4IGwgMC40NDE5Nzg5MjksMzYuNzgwNDY4OCBsIDAuMjg4NDYzMzA0LDM2LjE2NjQwNjMgbCAwLjI4ODQ2MzMwNCwzNS4yNTIxNDMzIGwgMC4yODg0NjMzMDQsMzQuMzk3MDg3NSBsIDAuMjg4NDYzMzA0LDMzLjYzMzM5ODQgYyAwLjI4ODQ2MzMwNCwzMy41MDA4NzA1IDAuNDQxOTc4OTI5LDMyLjk0MjU3ODEgMC40NDE5Nzg5MjksMzIuOTQyNTc4MSBsIDAuOTAyNTI1ODA0LDMyLjU1ODc4OTEgbCAxLjEzMjc5OTI0LDMxLjk0NDcyNjYgbCAxLjEzMjc5OTI0LDMxLjEwMDM5MDYgbCAxLjI4NjMxNDg3LDMwLjQwOTU3MDMgbCAxLjc0Njg2MTc0LDMwLjAyNTc4MTMgbCAxLjk3NzEzNTE4LDI5LjQxMTcxODggbCAyLjEzMDY1MDgsMjguNzIwODk4NCBsIDIuNTkxMTk3NjgsMjguMzM3MTA5NCBsIDIuOTc0OTg2NzQsMjcuODc2NTYyNSBsIDMuNDM1NTMzNjIsMjcuNDkyNzczNCBsIDMuODE5MzIyNjgsMjcuMDMyMjI2NiBsIDQuMjc5ODY5NTUsMjYuNjQ4NDM3NSBsIDQuNjYzNjU4NjIsMjYuMTg3ODkwNiBsIDUuMjc3NzIxMTIsMjUuOTU3NjE3MiBsIDYuMTIyMDU3MDUsMjUuOTU3NjE3MiBsIDYuODEyODc3MzcsMjUuODA0MTAxNiBsIDcuMTk2NjY2NDMsMjUuMzQzNTU0NyBsIDcuODEwNzI4OTMsMjUuMTEzMjgxMyBsIDguNjU1MDY0ODcsMjUuMTEzMjgxMyBsIDkuMzQ1ODg1MTgsMjQuOTU5NzY1NiBsIDkuNzI5Njc0MjQsMjQuNDk5MjE4OCBsIDEwLjM0MzczNjcsMjQuMjY4OTQ1MyBsIDEwLjkzNzA1ODUsMjQuMjY4OTQ1MyBsIDExLjY2OTUyOTQsMjQuMjY4OTQ1MyBsIDEyLjEyNDQ0MDQsMjQuMjY4OTQ1MyBsIDEyLjc4OTA5MjUsMjQuMjY4OTQ1MyBsIDEzLjcyMTA4MDUsMjQuMjY4OTQ1MyBjIDEzLjk1NjYwODQsMjQuMjY4OTQ1MyAxNC40MTE5MDA4LDI0LjExNTQyOTcgMTQuNDExOTAwOCwyNC4xMTU0Mjk3IGwgMTQuNzk1Njg5OSwyMy42NTQ4ODI4IGwgMTUuNDA5NzUyNCwyMy40MjQ2MDk0IGwgMTUuODkxNDkyOSwyMy40MjQ2MDk0IGwgMTYuMzc4MDM0MiwyMy40MjQ2MDk0IGwgMTYuNzE2NTA2OCwyMy40MjQ2MDk0IGwgMTcuMzI2NDYxOCwyMy40MjQ2MDk0IGwgMTcuOTQyNzYwMiwyMy40MjQ2MDk0IGMgMTguMTc4Mjg4MSwyMy40MjQ2MDk0IDE4LjYzMzU4MDUsMjMuMjcxMDkzOCAxOC42MzM1ODA1LDIzLjI3MTA5MzggbCAxOS4wMTczNjk2LDIyLjgxMDU0NjkgbCAxOS42MzE0MzIxLDIyLjU4MDI3MzQgbCAyMC4xMjU4ODc1LDIyLjU4MDI3MzQgbCAyMC42NTU2NjE2LDIyLjU4MDI3MzQgbCAyMS4zMDI2MzM4LDIyLjU4MDI3MzQgbCAyMi40OTE2NzMyLDIyLjU4MDI3MzQgbCAyMy40MzQ3ODMxLDIyLjU4MDI3MzQgbCAyMy44NTMxMTE3LDIyLjU4MDI3MzQgbCAyNC41NDM5MzIxLDIyLjQyNjc1NzggbCAyNC45Mjc3MjExLDIxLjk2NjIxMDkgbCAyNS4zODgyNjgsMjEuNTgyNDIxOSBsIDI1Ljc3MjA1NzEsMjEuMTIxODc1IGwgMjYuMzg2MTE5NiwyMC44OTE2MDE2IGwgMjcuMDc2OTM5OSwyMC43MzgwODU5IGwgMjcuMzA3MjEzMywyMC4xMjQwMjM0IGwgMjcuMzA3MjEzMywxOS4zMjAxMDY3IGwgMjcuMzA3MjEzMywxOC40MzUzNTE2IGwgMjcuMDc2OTM5OSwxNy43NDQ1MzEzIGwgMjYuNjE2MzkzLDE3LjM2MDc0MjIgbCAyNi42MTYzOTMsMTYuOTAwMTk1MyBsIDI3LjA3NjkzOTksMTYuNTE2NDA2MyBsIDI3LjA3NjkzOTksMTYuMDU1ODU5NCBsIDI2LjYxNjM5MywxNS42NzIwNzAzIGwgMjYuNDYyODc3NCwxNS4wNTgwMDc4IGwgMjYuMjMyNjAzOSwxNC4zNjcxODc1IGwgMjUuNzcyMDU3MSwxMy45ODMzOTg0IGwgMjUuMzg4MjY4LDEzLjUyMjg1MTYgbCAyNC45Mjc3MjExLDEzLjEzOTA2MjUgbCAyNC41NDM5MzIxLDEyLjY3ODUxNTYgbCAyMy44NTMxMTE3LDEyLjQ0ODI0MjIgbCAyMy4wMDg3NzU4LDEyLjQ0ODI0MjIgbCAyMi4zOTQ3MTMzLDEyLjI5NDcyNjYgbCAyMi4wMTA5MjQyLDExLjgzNDE3OTcgbCAyMS41NTAzNzc0LDExLjgzNDE3OTcgbCAyMS4xNjY1ODgzLDEyLjI5NDcyNjYgbCAyMC40NzU3NjgsMTIuNDQ4MjQyMiBsIDE5LjUzNTYxNTEsMTIuNDQ4MjQyMiBsIDE4LjY3ODkwMTgsMTIuNDQ4MjQyMiBsIDE4LjEwMDY0NjEsMTIuNDQ4MjQyMiBsIDE3LjQ5MDY5MTEsMTIuNDQ4MjQyMiBsIDE2Ljg2OTIwMjcsMTIuNDQ4MjQyMiBsIDE2LjI5NjI2NDksMTIuNDQ4MjQyMiBsIDE1Ljc1NTg1NTIsMTIuNDQ4MjQyMiBsIDE1LjE2NDQwODgsMTIuNDQ4MjQyMiBsIDE0LjU2NTQxNjQsMTIuNDQ4MjQyMiBsIDEzLjk1MTM1MzksMTIuNjc4NTE1NiBsIDEzLjU2NzU2NDksMTMuMTM5MDYyNSBsIDEzLjEwNzAxOCwxMy41MjI4NTE2IGwgMTIuNzIzMjI4OSwxMy45ODMzOTg0IGwgMTIuMjYyNjgyMSwxNC4zNjcxODc1IGwgMTIuMTA5MTY2NCwxNS4wNTgwMDc4IGwgMTIuMTA5MTY2NCwxNS43MDIwODk5IGwgMTIuMTA5MTY2NCwxNi42OTQ1NzkgbCAxMi4xMDkxNjY0LDE3LjU5MTAxNTYgbCAxMS44Nzg4OTMsMTguMjA1MDc4MSBsIDExLjQxODM0NjEsMTguNTg4ODY3MiBsIDExLjAzNDU1NzEsMTkuMDQ5NDE0MSBsIDEwLjU3NDAxMDIsMTkuNDMzMjAzMSBsIDEwLjE5MDIyMTEsMTkuODkzNzUgbCA5LjQ5OTQwMDgsMjAuMDQ3MjY1NiBsIDguNjQ0NDc4NCwyMC4wNDcyNjU2IGwgNy44MTA3Mjg5MywyMC4wNDcyNjU2IGwgNy4xOTY2NjY0MywxOS44OTM3NSBsIDYuODEyODc3MzcsMTkuNDMzMjAzMSBsIDYuMTIyMDU3MDUsMTkuMjAyOTI5NyBsIDUuNTA3OTk0NTUsMTkuMDQ5NDE0MSBsIDUuMTI0MjA1NDksMTguNTg4ODY3MiBsIDQuNjYzNjU4NjIsMTguMjA1MDc4MSBsIDQuNTEwMTQyOTksMTcuNTkxMDE1NiBsIDQuNTEwMTQyOTksMTYuNzQ2Njc5NyBsIDQuMjc5ODY5NTUsMTYuMDU1ODU5NCBsIDMuODE5MzIyNjgsMTUuNjcyMDcwMyBsIDMuNjY1ODA3MDUsMTUuMDU4MDA3OCBsIDMuNjY1ODA3MDUsMTQuMjEzNjcxOSBsIDMuODE5MzIyNjgsMTMuNTIyODUxNiBsIDQuMjc5ODY5NTUsMTMuMTM5MDYyNSBsIDQuNTEwMTQyOTksMTIuNTI1IGwgNC42NjM2NTg2MiwxMS44MzQxNzk3IGwgNS4xMjQyMDU0OSwxMS40NTAzOTA2IGwgNS41MDc5OTQ1NSwxMC45ODk4NDM4IGwgNS45Njg1NDE0MywxMC42MDYwNTQ3IGwgNi4zNTIzMzA0OSwxMC4xNDU1MDc4IGwgNi45NjYzOTI5OSw5LjkxNTIzNDM4IGwgNy42NTcyMTMzLDkuNzYxNzE4NzUgbCA4LjA0MTAwMjM3LDkuMzAxMTcxODggbCA4LjY1NTA2NDg3LDkuMDcwODk4NDQgbCA5LjM0NTg4NTE4LDguOTE3MzgyODEgbCA5LjcyOTY3NDI0LDguNDU2ODM1OTQgbCAxMC4zNDM3MzY3LDguMjI2NTYyNSBsIDEwLjk4OTA2MjEsOC4yMjY1NjI1IGwgMTEuNDM2MTAwMSw4LjIyNjU2MjUgbCAxMS45NDI4MDc1LDguMjI2NTYyNSBsIDEyLjI2ODE1ODMsOC4yMjY1NjI1IGwgMTIuOTU0OTc5Myw4LjIyNjU2MjUgbCAxMy43MjEwODA1LDguMjI2NTYyNSBsIDE0LjQxMTkwMDgsOC4wNzMwNDY4OCBsIDE0Ljc5NTY4OTksNy42MTI1IGwgMTUuMjU2MjM2Nyw3LjYxMjUgbCAxNS42NDAwMjU4LDguMDczMDQ2ODggbCAxNi4yNTQwODgzLDguMjI2NTYyNSBsIDE2Ljk2MDkxNyw4LjIyNjU2MjUgbCAxNy40OTE1MTk5LDguMjI2NTYyNSBsIDE4LjAyOTA5OCw4LjIyNjU2MjUgbCAxOC43OTMyNjgzLDguMjI2NTYyNSBsIDE5LjYzMTQzMjEsOC4yMjY1NjI1IGwgMjAuMzIyMjUyNCw4LjQ1NjgzNTk0IGwgMjAuNzA2MDQxNCw4LjkxNzM4MjgxIGwgMjEuMzIwMTAzOSw5LjA3MDg5ODQ0IGwgMjIuMDQ0MDgyOCw5LjA3MDg5ODQ0IGwgMjIuODc5NDU2LDkuMDcwODk4NDQgbCAyMy44NTMxMTE3LDkuMDcwODk4NDQgYyAyNC4wOTUwNDc0LDkuMDcwODk4NDQgMjQuNTQzOTMyMSw5LjMwMTE3MTg4IDI0LjU0MzkzMjEsOS4zMDExNzE4OCBsIDI0LjkyNzcyMTEsOS43NjE3MTg3NSBsIDI1LjU0MTc4MzYsOS45MTUyMzQzOCBsIDI2LjIzMjYwMzksMTAuMTQ1NTA3OCBsIDI2LjYxNjM5MywxMC42MDYwNTQ3IGwgMjcuMDc2OTM5OSwxMC45ODk4NDM4IGwgMjcuNDYwNzI4OSwxMS40NTAzOTA2IGwgMjcuOTIxMjc1OCwxMS44MzQxNzk3IGwgMjguMzA1MDY0OSwxMi4yOTQ3MjY2IGwgMjguNzY1NjExNywxMi42Nzg1MTU2IGwgMjkuMTQ5NDAwOCwxMy4xMzkwNjI1IGwgMjkuNjA5OTQ3NywxMy41MjI4NTE2IGwgMjkuODQwMjIxMSwxNC4yMTM2NzE5IGwgMjkuOTkzNzM2NywxNC44Mjc3MzQ0IGwgMzAuNDU0MjgzNiwxNS4yMTE1MjM0IGwgMzAuNjg0NTU3MSwxNS45MDIzNDM4IGwgMzAuNjg0NTU3MSwxNi41MDMyNzc0IGMgMzAuNjg0NTU3MSwxNi43MjcxNTMyIDMwLjY4NDU1NzEsMTYuOTUxMDI5IDMwLjY4NDU1NzEsMTcuMTc0OTA0OCBsIDMwLjY4NDU1NzEsMTcuODc2NTc0MSBsIDMwLjY4NDU1NzEsMTguNTQyMDU0OSBsIDMwLjY4NDU1NzEsMTkuMTY1MjAwOCBjIDMwLjY4NDU1NzEsMTkuNzA5MTE0NyAzMC42ODQ1NTcxLDIwLjEyNDAyMzQgMzAuNjg0NTU3MSwyMC4xMjQwMjM0IGwgMzAuODM4MDcyNywyMC43MzgwODU5IGwgMzEuMjk4NjE5NiwyMS4xMjE4NzUgbCAzMS41Mjg4OTMsMjEuODEyNjk1MyBsIDMxLjUyODg5MywyMi42MDYwNzM3IGMgMzEuNTI4ODkzLDIyLjg2Mjg2OSAzMS41Mjg4OTMsMjMuMTE5NjY0MyAzMS41Mjg4OTMsMjMuMzc2NDU5NiBsIDMxLjUyODg5MywyNC4xMTA1ODggYyAzMS41Mjg4OTMsMjQuMzc3MDc0MiAzMS41Mjg4OTMsMjQuNjU1NTg0NiAzMS41Mjg4OTMsMjQuOTQzNDc1IGwgMzEuNTI4ODkzLDI1Ljg4NzQxMzYgbCAzMS41Mjg4OTMsMjYuNzc0MTY4OSBjIDMxLjUyODg5MywyNy4wNzQ0MDQxIDMxLjUyODg5MywyNy4zNzQ2Mzk0IDMxLjUyODg5MywyNy42NzQ4NzQ3IGwgMzEuNTI4ODkzLDI4LjU4NzExMzggbCAzMS41Mjg4OTMsMjkuMjYyMTI1MiBsIDMxLjUyODg5MywzMC4wMTkyNTEyIGMgMzEuNTI4ODkzLDMwLjI2NDQ3NjggMzEuNTI4ODkzLDMwLjUxMjgwNjMgMzEuNTI4ODkzLDMwLjc2MTI1MjcgbCAzMS41Mjg4OTMsMzEuNTkxNjUzNCBsIDMxLjUyODg5MywzMi4yNzU2NDI5IGwgMzEuNTI4ODkzLDMzLjE1MDAzNjEgbCAzMS41Mjg4OTMsMzMuODEzODU5NCBsIDMxLjUyODg5MywzNC41NTMzNzQ2IGMgMzEuNTI4ODkzLDM1LjAzMTA3NzggMzEuNTI4ODkzLDM1LjMyMjA3MDMgMzEuNTI4ODkzLDM1LjMyMjA3MDMgbCAzMS42ODI0MDg2LDM1LjkzNjEzMjggbCAzMi4xNDI5NTU1LDM2LjMxOTkyMTkgbCAzMi41MjY3NDQ2LDM2Ljc4MDQ2ODggbCAzMi45ODcyOTE0LDM3LjE2NDI1NzggbCAzMy4zNzEwODA1LDM3LjYyNDgwNDcgbCAzMy45ODUxNDMsMzcuNzc4MzIwMyBsIDM0LjY3NTk2MzMsMzguMDA4NTkzOCBsIDM1LjA1OTc1MjQsMzguNDY5MTQwNiBsIDM1LjY3MzgxNDksMzguNjIyNjU2MyBsIDM2LjM2NDYzNTIsMzguODUyOTI5NyBsIDM2Ljc0ODQyNDIsMzkuMzEzNDc2NiBsIDM3LjM2MjQ4NjcsMzkuNDY2OTkyMiBsIDM4LjA1MzMwNzEsMzkuNjk3MjY1NiBsIDM4LjI4MzU4MDUsNDAuMzg4MDg1OSBsIDM4LjI4MzU4MDUsNDEuMDIxMjMyMyBsIDM4LjI4MzU4MDUsNDIuMDc1MzkzNyBsIDM4LjI4MzU4MDUsNDIuOTIxMDkzOCBsIDM4LjA1MzMwNzEsNDMuNTM1MTU2MyBsIDM3LjU5Mjc2MDIsNDMuOTE4OTQ1MyB6IG0gODcuNDE4MTg1Nyw0MC41NDE2MDE2IGwgODcuMjY0NjcwMSw0MS4yMzI0MjE5IGwgODcuMDM0Mzk2Niw0MS44NDY0ODQ0IGwgODYuNTczODQ5Nyw0Mi4yMzAyNzM0IGwgODYuMTkwMDYwNyw0Mi42OTA4MjAzIGwgODUuNzI5NTEzOCw0My4wNzQ2MDk0IGwgODUuMzQ1NzI0Nyw0My41MzUxNTYzIGwgODQuNjU0OTA0NCw0My42ODg2NzE5IGwgODQuMDQwODQxOSw0My45MTg5NDUzIGwgODMuNjU3MDUyOSw0NC4zNzk0OTIyIGwgODIuOTY2MjMyNiw0NC41MzMwMDc4IGwgODIuMzgyODE5MSw0NC41MzMwMDc4IGwgODEuNzYxNTQ2Myw0NC41MzMwMDc4IGwgODEuMTQ5Mzk5Myw0NC41MzMwMDc4IGwgODAuMzk5MTEyMiw0NC41MzMwMDc4IGwgNzkuNzY5ODc2Nyw0NC41MzMwMDc4IGwgNzkuMjY3NDE4OCw0NC41MzMwMDc4IGwgNzguNTk5MjY0Miw0NC41MzMwMDc4IGwgNzcuOTQ3MTI0Niw0NC41MzMwMDc4IGwgNzcuMjU5NzM0Miw0NC41MzMwMDc4IGwgNzYuNjY5NDE3Nyw0NC41MzMwMDc4IGwgNzYuMTU4OTk3LDQ0LjUzMzAwNzggbCA3NS41ODQ2MDYsNDQuNTMzMDA3OCBsIDc1LjAzNDE5MjYsNDQuNTMzMDA3OCBsIDc0LjUwNzc1Nyw0NC41MzMwMDc4IGwgNzQuMjcxMjAwOSw0NC41MzMwMDc4IGwgNzMuNzA5NjAzOSw0NC41MzMwMDc4IGwgNzIuODM0MjAxMyw0NC41MzMwMDc4IGMgNzIuODY4NjQxLDQ0LjUwMDA0MTMgNzIuMjIwMTM4OCw0NC4zNzk0OTIyIDcyLjIyMDEzODgsNDQuMzc5NDkyMiBsIDcxLjgzNjM0OTcsNDMuOTE4OTQ1MyBsIDcxLjM3NTgwMjksNDMuNTM1MTU2MyBsIDcxLjIyMjI4NzIsNDIuOTIxMDkzOCBsIDcwLjk5MjAxMzgsNDIuMjMwMjczNCBsIDcwLjUzMTQ2NjksNDEuODQ2NDg0NCBsIDcwLjM3Nzk1MTMsNDEuMjMyNDIxOSBsIDcwLjUzMTQ2NjksNDAuNTQxNjAxNiBsIDcwLjk5MjAxMzgsNDAuMTU3ODEyNSBsIDcxLjM3NTgwMjksMzkuNjk3MjY1NiBsIDcxLjgzNjM0OTcsMzkuMzEzNDc2NiBsIDcyLjIyMDEzODgsMzguODUyOTI5NyBsIDcyLjY4MDY4NTcsMzguNDY5MTQwNiBsIDczLjA2NDQ3NDcsMzguMDA4NTkzOCBsIDczLjY3ODUzNzIsMzcuNzc4MzIwMyBsIDc0LjM2OTM1NzYsMzcuNjI0ODA0NyBsIDc0Ljc1MzE0NjYsMzcuMTY0MjU3OCBsIDc1LjIxMzY5MzUsMzYuNzgwNDY4OCBsIDc1LjU5NzQ4MjYsMzYuMzE5OTIxOSBsIDc2LjA1ODAyOTQsMzUuOTM2MTMyOCBsIDc2LjI4ODMwMjksMzUuMzIyMDcwMyBsIDc2LjI4ODMwMjksMzQuNDc3NzM0NCBsIDc2LjQ0MTgxODUsMzMuNzg2OTE0MSBsIDc2LjkwMjM2NTQsMzMuNDAzMTI1IGwgNzcuMTMyNjM4OCwzMi43ODkwNjI1IGwgNzcuMTMyNjM4OCwzMi4wNTg2NTE3IGwgNzcuMTMyNjM4OCwzMS4xMzQ3OTQ3IGwgNzcuMTMyNjM4OCwzMC4wNjgwNTU3IGwgNzcuMTMyNjM4OCwyOS4yNjMwMTM2IGwgNzcuMTMyNjM4OCwyOC4zODcwMjI1IGwgNzcuMTMyNjM4OCwyNy41MDQ5NDc1IGwgNzcuMTMyNjM4OCwyNi42MDI4MzE1IGwgNzcuMTMyNjM4OCwyNS43NzE2NjQ0IGwgNzcuMTMyNjM4OCwyNC44MjM3NDAzIGwgNzcuMTMyNjM4OCwyMy45MjE3MTM4IGwgNzcuMTMyNjM4OCwyMy4wOTA0NTcyIGwgNzcuMTMyNjM4OCwyMi4xMzEzNDk0IGwgNzcuMTMyNjM4OCwyMS4xNjQyNzg5IGwgNzcuMTMyNjM4OCwyMC4yMzgyNzQ3IGwgNzcuMTMyNjM4OCwxOS41NDI5MjE1IGwgNzcuMTMyNjM4OCwxOC45MjE2NDg4IGwgNzcuMTMyNjM4OCwxOC4xNzk0MTM5IGwgNzcuMTMyNjM4OCwxNy41OTEwMTU2IGwgNzcuMjg2MTU0NCwxNi45MDAxOTUzIGwgNzcuNzQ2NzAxMywxNi41MTY0MDYzIGwgNzcuOTc2OTc0NywxNS45MDIzNDM4IGwgNzcuNzQ2NzAxMywxNS4yMTE1MjM0IGwgNzcuMjg2MTU0NCwxNC44Mjc3MzQ0IGwgNzcuMTMyNjM4OCwxNC4yMTM2NzE5IGwgNzYuOTAyMzY1NCwxMy41MjI4NTE2IGwgNzYuNDQxODE4NSwxMy4xMzkwNjI1IGwgNzYuMDU4MDI5NCwxMi42Nzg1MTU2IGwgNzUuMzY3MjA5MSwxMi40NDgyNDIyIGwgNzQuNzUzMTQ2NiwxMi4yOTQ3MjY2IGwgNzQuMzY5MzU3NiwxMS44MzQxNzk3IGwgNzMuNjc4NTM3MiwxMS42MDM5MDYzIGwgNzMuMDY0NDc0NywxMS40NTAzOTA2IGwgNzIuNjgwNjg1NywxMC45ODk4NDM4IGwgNzEuOTg5ODY1NCwxMC43NTk1NzAzIGwgNzEuMjEzNTA4MSwxMC43NTk1NzAzIGwgNzAuMjg4NDg4LDEwLjc1OTU3MDMgbCA2OS43NzExNzgyLDEwLjc1OTU3MDMgbCA2OS4wNTM5OTQ2LDEwLjc1OTU3MDMgbCA2OC4yOTY4MTg0LDEwLjc1OTU3MDMgbCA2Ny41NjU4NTY2LDEwLjc1OTU3MDMgbCA2Ni45MDQ1OTEyLDEwLjc1OTU3MDMgbCA2Ni4wNzk1MTM4LDEwLjc1OTU3MDMgYyA2Ni4xMTYzODkzLDEwLjgxMDIwMSA2NS40NjU0NTEzLDEwLjk4OTg0MzggNjUuNDY1NDUxMywxMC45ODk4NDM4IGwgNjUuMDgxNjYyMiwxMS40NTAzOTA2IGwgNjQuMzkwODQxOSwxMS42MDM5MDYzIGwgNjMuNzc2Nzc5NCwxMS44MzQxNzk3IGwgNjMuNjIzMjYzOCwxMi41MjUgbCA2My4zOTI5OTA0LDEzLjEzOTA2MjUgbCA2Mi45MzI0NDM1LDEzLjUyMjg1MTYgbCA2Mi41NDg2NTQ0LDEzLjk4MzM5ODQgbCA2Mi4wODgxMDc2LDE0LjM2NzE4NzUgbCA2MS43MDQzMTg1LDE0LjgyNzczNDQgbCA2MS4yNDM3NzE2LDE1LjIxMTUyMzQgbCA2MS4wOTAyNTYsMTUuOTAyMzQzOCBsIDYwLjg1OTk4MjYsMTYuNTE2NDA2MyBsIDYwLjM5OTQzNTcsMTYuOTAwMTk1MyBsIDYwLjAxNTY0NjYsMTcuMzYwNzQyMiBsIDU5LjU1NTA5OTcsMTcuNzQ0NTMxMyBsIDU5LjQwMTU4NDEsMTguNDM1MzUxNiBsIDU5LjE3MTMxMDcsMTkuMDQ5NDE0MSBsIDU4LjcxMDc2MzgsMTkuNDMzMjAzMSBsIDU4LjU1NzI0ODIsMjAuMTI0MDIzNCBsIDU4LjMyNjk3NDcsMjAuNzM4MDg1OSBsIDU3Ljg2NjQyNzksMjEuMTIxODc1IGwgNTcuNzEyOTEyMiwyMS44MTI2OTUzIGwgNTcuNDgyNjM4OCwyMi40MjY3NTc4IGwgNTcuMDIyMDkxOSwyMi44MTA1NDY5IGwgNTYuODY4NTc2MywyMy41MDEzNjcyIGwgNTYuODY4NTc2MywyMy45ODk3MTAyIGwgNTYuODY4NTc2MywyNC42OTg5MzEgbCA1Ni44Njg1NzYzLDI1LjY1NDgxNzkgbCA1Ni44Njg1NzYzLDI2LjMyNTExOTcgbCA1Ni44Njg1NzYzLDI3LjA5MDI1ODcgbCA1Ni44Njg1NzYzLDI3LjgwMDQ2MzcgbCA1Ni44Njg1NzYzLDI4LjU2NzM4MjggbCA1Ny4wMjIwOTE5LDI5LjE4MTQ0NTMgbCA1Ny40ODI2Mzg4LDI5LjU2NTIzNDQgbCA1Ny43MTI5MTIyLDMwLjI1NjA1NDcgbCA1Ny43MTI5MTIyLDMxLjEwMDM5MDYgbCA1Ny40ODI2Mzg4LDMxLjcxNDQ1MzEgbCA1Ny4wMjIwOTE5LDMyLjA5ODI0MjIgbCA1Ni44Njg1NzYzLDMyLjc4OTA2MjUgbCA1Ny4wMjIwOTE5LDMzLjQwMzEyNSBsIDU3LjQ4MjYzODgsMzMuNzg2OTE0MSBsIDU3LjcxMjkxMjIsMzQuNDc3NzM0NCBsIDU3LjcxMjkxMjIsMzUuMzIyMDcwMyBsIDU3Ljg2NjQyNzksMzUuOTM2MTMyOCBsIDU4LjMyNjk3NDcsMzYuMzE5OTIxOSBsIDU4LjcxMDc2MzgsMzYuNzgwNDY4OCBsIDU5LjE3MTMxMDcsMzcuMTY0MjU3OCBsIDU5LjU1NTA5OTcsMzcuNjI0ODA0NyBsIDYwLjE2OTE2MjIsMzcuNzc4MzIwMyBsIDYwLjg1OTk4MjYsMzguMDA4NTkzOCBsIDYxLjI0Mzc3MTYsMzguNDY5MTQwNiBsIDYxLjcwNDMxODUsMzguODUyOTI5NyBsIDYyLjA4ODEwNzYsMzkuMzEzNDc2NiBsIDYyLjU0ODY1NDQsMzkuNjk3MjY1NiBsIDYyLjc3ODkyNzksNDAuMzg4MDg1OSBsIDYyLjkzMjQ0MzUsNDEuMDAyMTQ4NCBsIDYzLjM5Mjk5MDQsNDEuMzg1OTM3NSBsIDYzLjM5Mjk5MDQsNDEuODQ2NDg0NCBsIDYyLjkzMjQ0MzUsNDIuMjMwMjczNCBsIDYyLjc3ODkyNzksNDIuOTIxMDkzOCBsIDYyLjU0ODY1NDQsNDMuNTM1MTU2MyBsIDYxLjg1NzgzNDEsNDMuNjg4NjcxOSBsIDYxLjI0Mzc3MTYsNDMuOTE4OTQ1MyBsIDYwLjg1OTk4MjYsNDQuMzc5NDkyMiBsIDYwLjE2OTE2MjIsNDQuNTMzMDA3OCBsIDU5LjMyNDgyNjMsNDQuNTMzMDA3OCBsIDU4LjcxMDc2MzgsNDQuNzYzMjgxMyBsIDU4LjMyNjk3NDcsNDUuMjIzODI4MSBsIDU3Ljg2NjQyNzksNDUuMjIzODI4MSBsIDU3LjQ4MjYzODgsNDQuNzYzMjgxMyBsIDU3LjAyMjA5MTksNDQuNzYzMjgxMyBsIDU2LjYzODMwMjksNDUuMjIzODI4MSBsIDU1Ljk0NzQ4MjYsNDUuMzc3MzQzOCBsIDU1LjEwMzE0NjYsNDUuMzc3MzQzOCBsIDU0LjQ4OTA4NDEsNDUuMjIzODI4MSBsIDU0LjEwNTI5NTEsNDQuNzYzMjgxMyBsIDUzLjY0NDc0ODIsNDQuNzYzMjgxMyBsIDUzLjI2MDk1OTEsNDUuMjIzODI4MSBsIDUyLjU3MDEzODgsNDUuMzc3MzQzOCBsIDUxLjQ5NzQ5MTEsNDUuMzc3MzQzOCBsIDUwLjU5MDYzMzIsNDUuMzc3MzQzOCBsIDQ5Ljk4NTM3NTQsNDUuMzc3MzQzOCBsIDQ5LjE5Mjc5NTEsNDUuMzc3MzQzOCBjIDQ5LjIyNjg0MDIsNDUuMzQ0MzAzOSA0OC41Nzg3MzI2LDQ1LjIyMzgyODEgNDguNTc4NzMyNiw0NS4yMjM4MjgxIGwgNDguMTk0OTQzNSw0NC43NjMyODEzIGwgNDcuNTA0MTIzMiw0NC41MzMwMDc4IGwgNDYuODkwMDYwNyw0NC4zNzk0OTIyIGwgNDYuNTA2MjcxNiw0My45MTg5NDUzIGwgNDYuMDQ1NzI0Nyw0My41MzUxNTYzIGwgNDUuNjYxOTM1Nyw0My4wNzQ2MDk0IGwgNDUuMjAxMzg4OCw0Mi42OTA4MjAzIGwgNDUuMDQ3ODczMiw0Mi4wNzY3NTc4IGwgNDUuMDQ3ODczMiw0MS4yMzI0MjE5IGwgNDUuMjAxMzg4OCw0MC41NDE2MDE2IGwgNDUuNjYxOTM1Nyw0MC4xNTc4MTI1IGwgNDYuMDQ1NzI0NywzOS42OTcyNjU2IGwgNDYuNTA2MjcxNiwzOS4zMTM0NzY2IGwgNDYuODkwMDYwNywzOC44NTI5Mjk3IGwgNDcuNTA0MTIzMiwzOC42MjI2NTYzIGwgNDguMTk0OTQzNSwzOC40NjkxNDA2IGwgNDguNTc4NzMyNiwzOC4wMDg1OTM4IGwgNDkuMDM5Mjc5NCwzNy42MjQ4MDQ3IGwgNDkuNDIzMDY4NSwzNy4xNjQyNTc4IGwgNDkuODgzNjE1NCwzNi43ODA0Njg4IGwgNTAuMjY3NDA0NCwzNi4zMTk5MjE5IGwgNTAuNzI3OTUxMywzNS45MzYxMzI4IGwgNTAuOTU4MjI0NywzNS4zMjIwNzAzIGwgNTAuOTU4MjI0NywzNC42OTk4NjYzIGwgNTAuOTU4MjI0NywzMy45MTUwNDQxIGwgNTAuOTU4MjI0NywzMi45NDAwMTA5IGwgNTAuOTU4MjI0NywzMi4zMDc1NTQ0IGwgNTAuOTU4MjI0NywzMS40NjgzMzUyIGwgNTAuOTU4MjI0NywzMC41ODkxMjMyIGwgNTAuOTU4MjI0NywyOS44MTM4NzQyIGwgNTAuOTU4MjI0NywyOS4wMjM2ODQgbCA1MC45NTgyMjQ3LDI4LjEyMzgwNDcgbCA1MC45NTgyMjQ3LDI3LjAwMzkyMTEgbCA1MC45NTgyMjQ3LDI2LjAzNDM3NSBjIDUwLjk1ODIyNDksMjUuODQ2MDEyNSA1MS4xMTE3NDA0LDI1LjM0MzU1NDcgNTEuMTExNzQwNCwyNS4zNDM1NTQ3IGwgNTEuNTcyMjg3MiwyNC45NTk3NjU2IGwgNTEuODAyNTYwNywyNC4zNDU3MDMxIGwgNTEuODAyNTYwNywyMy41MDEzNjcyIGwgNTEuNTcyMjg3MiwyMi44MTA1NDY5IGwgNTEuMTExNzQwNCwyMi40MjY3NTc4IGwgNTAuOTU4MjI0NywyMS44MTI2OTUzIGwgNTAuOTU4MjI0NywyMS4wMjkwOTEyIGwgNTAuOTU4MjI0NywyMC4zMzU4ODU0IGwgNTAuOTU4MjI0NywxOS42MzQ3MTY3IGwgNTAuOTU4MjI0NywxOC45NTc1MjU4IGwgNTAuOTU4MjI0NywxOC4wMDA1NjUzIGwgNTAuOTU4MjI0NywxNy4zNTUzMTQ4IGwgNTAuOTU4MjI0NywxNi43NDY2Nzk3IGwgNTAuNzI3OTUxMywxNi4wNTU4NTk0IGwgNTAuMjY3NDA0NCwxNS42NzIwNzAzIGwgNTAuMTEzODg4OCwxNS4wNTgwMDc4IGwgNDkuODgzNjE1NCwxNC4zNjcxODc1IGwgNDkuNDIzMDY4NSwxMy45ODMzOTg0IGwgNDkuMDM5Mjc5NCwxMy41MjI4NTE2IGwgNDguNTc4NzMyNiwxMy4xMzkwNjI1IGwgNDguMTk0OTQzNSwxMi42Nzg1MTU2IGwgNDcuNzM0Mzk2NiwxMi4yOTQ3MjY2IGwgNDcuMzUwNjA3NiwxMS44MzQxNzk3IGwgNDYuODkwMDYwNywxMS40NTAzOTA2IGwgNDYuNTA2MjcxNiwxMC45ODk4NDM4IGwgNDUuODE1NDUxMywxMC43NTk1NzAzIGwgNDUuMjAxMzg4OCwxMC42MDYwNTQ3IGwgNDUuMjAxMzg4OCwxMC4xNDU1MDc4IGwgNDUuNjYxOTM1Nyw5Ljc2MTcxODc1IGwgNDUuODkyMjA5MSw5LjE0NzY1NjI1IGwgNDUuODkyMjA5MSw4LjMwMzMyMDMxIGwgNDYuMDQ1NzI0Nyw3LjYxMjUgbCA0Ni41MDYyNzE2LDcuMjI4NzEwOTQgbCA0Ni44OTAwNjA3LDYuNzY4MTY0MDYgbCA0Ny4zNTA2MDc2LDYuMzg0Mzc1IGwgNDcuNzM0Mzk2Niw1LjkyMzgyODEzIGwgNDguMzQ4NDU5MSw1LjY5MzU1NDY5IGwgNDkuMDM5Mjc5NCw1LjU0MDAzOTA2IGwgNDkuNDIzMDY4NSw1LjA3OTQ5MjE5IGwgNTAuMDM3MTMxLDQuODQ5MjE4NzUgbCA1MC43Mjc5NTEzLDUuMDc5NDkyMTkgbCA1MS4xMTE3NDA0LDUuNTQwMDM5MDYgbCA1MS43MjU4MDI5LDUuNjkzNTU0NjkgbCA1Mi41NzAxMzg4LDUuNjkzNTU0NjkgbCA1My4yNjA5NTkxLDUuOTIzODI4MTMgbCA1My42NDQ3NDgyLDYuMzg0Mzc1IGwgNTQuMjU4ODEwNyw2LjUzNzg5MDYzIGwgNTQuOTQ5NjMxLDYuNzY4MTY0MDYgbCA1NS4zMzM0MjAxLDcuMjI4NzEwOTQgbCA1NS43OTM5NjY5LDcuNjEyNSBsIDU2LjE3Nzc1Niw4LjA3MzA0Njg4IGwgNTYuNjM4MzAyOSw4LjQ1NjgzNTk0IGwgNTcuMDIyMDkxOSw4LjkxNzM4MjgxIGwgNTcuNDgyNjM4OCw5LjMwMTE3MTg4IGwgNTcuODY2NDI3OSw5Ljc2MTcxODc1IGwgNTguMzI2OTc0NywxMC4xNDU1MDc4IGwgNTguNzEwNzYzOCwxMC42MDYwNTQ3IGwgNTkuMTcxMzEwNywxMC42MDYwNTQ3IGwgNTkuNTU1MDk5NywxMC4xNDU1MDc4IGwgNjAuMTY5MTYyMiw5LjkxNTIzNDM4IGwgNjEuMDEzNDk4Miw5LjkxNTIzNDM4IGwgNjEuNzA0MzE4NSw5Ljc2MTcxODc1IGwgNjIuMDg4MTA3Niw5LjMwMTE3MTg4IGwgNjIuNTQ4NjU0NCw4LjkxNzM4MjgxIGwgNjIuOTMyNDQzNSw4LjQ1NjgzNTk0IGwgNjMuMzkyOTkwNCw4LjA3MzA0Njg4IGwgNjMuNzc2Nzc5NCw3LjYxMjUgbCA2NC4zOTA4NDE5LDcuMzgyMjI2NTYgbCA2NS4wODE2NjIyLDcuMjI4NzEwOTQgbCA2NS40NjU0NTEzLDYuNzY4MTY0MDYgbCA2Ni4wNzk1MTM4LDYuNTM3ODkwNjMgbCA2Ni42NDU2Njc5LDYuNTM3ODkwNjMgbCA2Ny4zMDIxMDIsNi41Mzc4OTA2MyBsIDY4LjA1OTE4ODcsNi41Mzc4OTA2MyBsIDY4LjYxMjUyMTYsNi41Mzc4OTA2MyBsIDY5LjMwMzM0MTksNi4zODQzNzUgbCA2OS42ODcxMzEsNS45MjM4MjgxMyBsIDcwLjMwMTE5MzUsNS42OTM1NTQ2OSBsIDcwLjk5MDgxOTcsNS42OTM1NTQ2OSBsIDcxLjU0ODEyMjIsNS42OTM1NTQ2OSBsIDcyLjA0MzYwMTUsNS42OTM1NTQ2OSBsIDcyLjc4NTgzNjQsNS42OTM1NTQ2OSBsIDczLjY3ODUzNzIsNS42OTM1NTQ2OSBsIDc0LjM2OTM1NzYsNS45MjM4MjgxMyBsIDc0Ljc1MzE0NjYsNi4zODQzNzUgbCA3NS4zNjcyMDkxLDYuNTM3ODkwNjMgbCA3Ni4wNTgwMjk0LDYuNzY4MTY0MDYgbCA3Ni40NDE4MTg1LDcuMjI4NzEwOTQgbCA3Ny4wNTU4ODEsNy4zODIyMjY1NiBsIDc3Ljc0NjcwMTMsNy42MTI1IGwgNzguMTMwNDkwNCw4LjA3MzA0Njg4IGwgNzguNTkxMDM3Miw4LjQ1NjgzNTk0IGwgNzguOTc0ODI2Myw4LjkxNzM4MjgxIGwgNzkuNDM1MzczMiw5LjMwMTE3MTg4IGwgNzkuODE5MTYyMiw5Ljc2MTcxODc1IGwgODAuMjc5NzA5MSwxMC4xNDU1MDc4IGwgODAuNTA5OTgyNiwxMC44MzYzMjgxIGwgODAuNjYzNDk4MiwxMS40NTAzOTA2IGwgODEuMTI0MDQ1MSwxMS44MzQxNzk3IGwgODEuMzU0MzE4NSwxMi41MjUgbCA4MS41MDc4MzQxLDEzLjEzOTA2MjUgbCA4MS45NjgzODEsMTMuNTIyODUxNiBsIDgyLjE5ODY1NDQsMTQuMjEzNjcxOSBsIDgyLjE5ODY1NDQsMTQuOTQwMDk4OSBsIDgyLjE5ODY1NDQsMTUuOTc5MTAyNSBsIDgyLjE5ODY1NDQsMTYuNzQ2Njc5NyBsIDgxLjk2ODM4MSwxNy4zNjA3NDIyIGwgODEuNTA3ODM0MSwxNy43NDQ1MzEzIGwgODEuNTA3ODM0MSwxOC4yMDUwNzgxIGwgODEuOTY4MzgxLDE4LjU4ODg2NzIgbCA4MS45NjgzODEsMTkuMDQ5NDE0MSBsIDgxLjUwNzgzNDEsMTkuNDMzMjAzMSBsIDgxLjM1NDMxODUsMjAuMTI0MDIzNCBsIDgxLjM1NDMxODUsMjAuOTY1Mzg5MyBsIDgxLjM1NDMxODUsMjEuODEyNjk1MyBjIDgxLjM1NDMxOTMsMjIuNDI2NzU4MiA4MS41MDc4MzQxLDIyLjQyNjc1NzggODEuNTA3ODM0MSwyMi40MjY3NTc4IGwgODEuOTY4MzgxLDIyLjgxMDU0NjkgbCA4Mi4xOTg2NTQ0LDIzLjUwMTM2NzIgbCA4Mi4xOTg2NTQ0LDI0LjM0MjMwNzYgbCA4Mi4xOTg2NTQ0LDI1LjIxOTI4MjggbCA4Mi4xOTg2NTQ0LDI1LjkyOTU3NzIgbCA4Mi4xOTg2NTQ0LDI3LjA3MTM4MDcgbCA4Mi4xOTg2NTQ0LDI3LjY4NDY5MDggbCA4Mi4xOTg2NTQ0LDI4LjU2NzM4MjggYyA4Mi4xOTg2NTI5LDI5LjE4MTQ0NDUgODIuMzUyMTcwMSwyOS4xODE0NDUzIDgyLjM1MjE3MDEsMjkuMTgxNDQ1MyBsIDgyLjgxMjcxNjksMjkuNTY1MjM0NCBsIDgzLjA0Mjk5MDQsMzAuMjU2MDU0NyBsIDgzLjA0Mjk5MDQsMzEuMTAwMzkwNiBsIDgyLjgxMjcxNjksMzEuNzE0NDUzMSBsIDgyLjM1MjE3MDEsMzIuMDk4MjQyMiBsIDgyLjE5ODY1NDQsMzIuNzg5MDYyNSBsIDgyLjM1MjE3MDEsMzMuNDAzMTI1IGwgODIuODEyNzE2OSwzMy43ODY5MTQxIGwgODMuMDQyOTkwNCwzNC40Nzc3MzQ0IGwgODMuMTk2NTA2LDM1LjA5MTc5NjkgbCA4My42NTcwNTI5LDM1LjQ3NTU4NTkgbCA4NC4wNDA4NDE5LDM1LjkzNjEzMjggbCA4NC42NTQ5MDQ0LDM2LjA4OTY0ODQgbCA4NS4zNDU3MjQ3LDM2LjMxOTkyMTkgbCA4NS43Mjk1MTM4LDM2Ljc4MDQ2ODggbCA4Ni4zNDM1NzYzLDM2LjkzMzk4NDQgbCA4Ny4wMzQzOTY2LDM3LjE2NDI1NzggbCA4Ny4yNjQ2NzAxLDM3Ljg1NTA3ODEgbCA4Ny40MTgxODU3LDM4LjQ2OTE0MDYgbCA4Ny44Nzg3MzI2LDM4Ljg1MjkyOTcgbCA4OC4xMDkwMDYsMzkuNTQzNzUgbCA4Ny44Nzg3MzI2LDQwLjE1NzgxMjUgbCA4Ny40MTgxODU3LDQwLjU0MTYwMTYgeiBtIDEzMC40ODg5MjQsMzIuMDk4MjQyMiBsIDEzMC4zMzU0MDgsMzIuNzg5MDYyNSBsIDEzMC4zMzU0MDgsMzMuNjg4NjYwNiBsIDEzMC4zMzU0MDgsMzQuNDc3NzM0NCBsIDEzMC4xMDUxMzUsMzUuMDkxNzk2OSBsIDEyOS42NDQ1ODgsMzUuNDc1NTg1OSBsIDEyOS40OTEwNzIsMzYuMTY2NDA2MyBsIDEyOS40OTEwNzIsMzYuOTc0MjYzMSBsIDEyOS40OTEwNzIsMzcuODU1MDc4MSBsIDEyOS4yNjA3OTksMzguNDY5MTQwNiBsIDEyOC44MDAyNTIsMzguODUyOTI5NyBsIDEyOC42NDY3MzYsMzkuNTQzNzUgbCAxMjguNDE2NDYzLDQwLjE1NzgxMjUgbCAxMjcuOTU1OTE2LDQwLjU0MTYwMTYgbCAxMjcuNTcyMTI3LDQxLjAwMjE0ODQgbCAxMjcuMTExNTgsNDEuMzg1OTM3NSBsIDEyNi43Mjc3OTEsNDEuODQ2NDg0NCBsIDEyNi4yNjcyNDQsNDIuMjMwMjczNCBsIDEyNS44ODM0NTUsNDIuNjkwODIwMyBsIDEyNS40MjI5MDgsNDMuMDc0NjA5NCBsIDEyNS4wMzkxMTksNDMuNTM1MTU2MyBsIDEyNC41Nzg1NzIsNDMuOTE4OTQ1MyBsIDEyNC4xOTQ3ODMsNDQuMzc5NDkyMiBsIDEyMy41MDM5NjMsNDQuNTMzMDA3OCBsIDEyMi44ODk5LDQ0Ljc2MzI4MTMgbCAxMjIuNTA2MTExLDQ1LjIyMzgyODEgbCAxMjEuODE1MjkxLDQ1LjM3NzM0MzggbCAxMjAuOTcwOTU1LDQ1LjM3NzM0MzggbCAxMjAuMzU2ODkyLDQ1LjYwNzYxNzIgbCAxMTkuOTczMTAzLDQ2LjA2ODE2NDEgbCAxMTkuMjgyMjgzLDQ2LjIyMTY3OTcgbCAxMTguNjI0MTg1LDQ2LjIyMTY3OTcgbCAxMTcuOTQxNTY2LDQ2LjIyMTY3OTcgbCAxMTcuMTE4NTA3LDQ2LjIyMTY3OTcgbCAxMTYuNDQyMDQzLDQ2LjIyMTY3OTcgbCAxMTUuNzA3NTQ4LDQ2LjIyMTY3OTcgbCAxMTUuMDYwNjAzLDQ2LjIyMTY3OTcgbCAxMTQuNDQ2NTQxLDQ2LjA2ODE2NDEgbCAxMTQuMDYyNzUyLDQ1LjYwNzYxNzIgbCAxMTMuMzcxOTMyLDQ1LjM3NzM0MzggbCAxMTIuNTI3NTk2LDQ1LjM3NzM0MzggbCAxMTEuOTEzNTMzLDQ1LjIyMzgyODEgbCAxMTEuNTI5NzQ0LDQ0Ljc2MzI4MTMgbCAxMTAuODM4OTI0LDQ0LjUzMzAwNzggbCAxMTAuMjI0ODYxLDQ0LjM3OTQ5MjIgbCAxMDkuODQxMDcyLDQzLjkxODk0NTMgbCAxMDkuMzgwNTI1LDQzLjUzNTE1NjMgbCAxMDguOTk2NzM2LDQzLjA3NDYwOTQgbCAxMDguNTM2MTg5LDQyLjY5MDgyMDMgbCAxMDguMTUyNCw0Mi4yMzAyNzM0IGwgMTA3LjY5MTg1Myw0MS44NDY0ODQ0IGwgMTA3LjMwODA2NCw0MS4zODU5Mzc1IGwgMTA2Ljg0NzUxNyw0MS4wMDIxNDg0IGwgMTA2LjQ2MzcyOCw0MC41NDE2MDE2IGwgMTA2LjAwMzE4Miw0MC4xNTc4MTI1IGwgMTA1LjYxOTM5MiwzOS42OTcyNjU2IGwgMTA1LjE1ODg0NiwzOS4zMTM0NzY2IGwgMTA1LjAwNTMzLDM4LjY5OTQxNDEgbCAxMDQuNzc1MDU3LDM4LjAwODU5MzggbCAxMDQuMzE0NTEsMzcuNjI0ODA0NyBsIDEwNC4xNjA5OTQsMzcuMDEwNzQyMiBsIDEwNC4xNjA5OTQsMzYuMTY2NDA2MyBsIDEwMy45MzA3MjEsMzUuNDc1NTg1OSBsIDEwMy40NzAxNzQsMzUuMDkxNzk2OSBsIDEwMy4zMTY2NTgsMzQuNDc3NzM0NCBsIDEwMy4zMTY2NTgsMzMuNzQ0MDUzNCBsIDEwMy4zMTY2NTgsMzIuNzg5MDYyNSBjIDEwMy4zMTY2NTgsMzIuNDMxMTcwNyAxMDMuNDcwMTc0LDMyLjA5ODI0MjIgMTAzLjQ3MDE3NCwzMi4wOTgyNDIyIGwgMTAzLjkzMDcyMSwzMS43MTQ0NTMxIGwgMTA0LjE2MDk5NCwzMS4xMDAzOTA2IGwgMTA0LjE2MDk5NCwzMC4zOTYyNjM5IGwgMTA0LjE2MDk5NCwyOS41MTUwOTQxIGwgMTA0LjE2MDk5NCwyOC40MzA2NTc0IGwgMTA0LjE2MDk5NCwyNy43MjMwNDY5IGMgMTA0LjE2MDk5NCwyNy41NDQyOTE3IDEwNC4zMTQ1MSwyNy4wMzIyMjY2IDEwNC4zMTQ1MSwyNy4wMzIyMjY2IGwgMTA0Ljc3NTA1NywyNi42NDg0Mzc1IGwgMTA1LjAwNTMzLDI2LjAzNDM3NSBsIDEwNC43NzUwNTcsMjUuMzQzNTU0NyBsIDEwNC4zMTQ1MSwyNC45NTk3NjU2IGwgMTA0LjE2MDk5NCwyNC4zNDU3MDMxIGwgMTA0LjMxNDUxLDIzLjY1NDg4MjggbCAxMDQuNzc1MDU3LDIzLjI3MTA5MzggbCAxMDUuMDA1MzMsMjIuNjU3MDMxMyBsIDEwNS4wMDUzMywyMS41NDExNjM0IGwgMTA1LjAwNTMzLDIwLjYyNzE0MTYgbCAxMDUuMDA1MzMsMTkuOTkyMjQyMiBsIDEwNS4wMDUzMywxOS4yNTY3MDg1IGwgMTA1LjAwNTMzLDE4LjQ5NTU5NjUgbCAxMDUuMDA1MzMsMTcuNjQ3NzU4NCBsIDEwNS4wMDUzMywxNi44MDAwMDAyIGwgMTA1LjAwNTMzLDE1LjkwMjM0MzggbCAxMDQuNzc1MDU3LDE1LjIxMTUyMzQgbCAxMDQuMzE0NTEsMTQuODI3NzM0NCBsIDEwNC4xNjA5OTQsMTQuMjEzNjcxOSBsIDEwMy45MzA3MjEsMTMuNTIyODUxNiBsIDEwMy40NzAxNzQsMTMuMTM5MDYyNSBsIDEwMy4wODYzODUsMTIuNjc4NTE1NiBsIDEwMi4zOTU1NjQsMTIuNDQ4MjQyMiBsIDEwMS43ODE1MDIsMTIuMjk0NzI2NiBsIDEwMS4zOTc3MTMsMTEuODM0MTc5NyBsIDEwMC43MDY4OTIsMTEuNjAzOTA2MyBsIDEwMC4wNTMzNTcsMTEuNjAzOTA2MyBsIDk5LjQyODg0ODMsMTEuNjAzOTA2MyBsIDk4LjgyMDAwNjcsMTEuNjAzOTA2MyBsIDk4LjE3Mzg4NDYsMTEuNjAzOTA2MyBsIDk3LjU1OTgyMjEsMTEuNDUwMzkwNiBsIDk3LjE3NjAzMzEsMTAuOTg5ODQzOCBsIDk2LjcxNTQ4NjIsMTAuNjA2MDU0NyBsIDk2LjMzMTY5NzEsMTAuMTQ1NTA3OCBsIDk1Ljg3MTE1MDMsOS43NjE3MTg3NSBsIDk1LjQ4NzM2MTIsOS4zMDExNzE4OCBsIDk1LjAyNjgxNDMsOC45MTczODI4MSBsIDk0Ljg3MzI5ODcsOC4zMDMzMjAzMSBsIDk1LjAyNjgxNDMsNy42MTI1IGwgOTUuNDg3MzYxMiw3LjIyODcxMDk0IGwgOTUuODcxMTUwMyw2Ljc2ODE2NDA2IGwgOTYuMzMxNjk3MSw2LjM4NDM3NSBsIDk2LjcxNTQ4NjIsNS45MjM4MjgxMyBsIDk3LjE3NjAzMzEsNS41NDAwMzkwNiBsIDk3LjU1OTgyMjEsNS4wNzk0OTIxOSBsIDk4LjE3Mzg4NDYsNC44NDkyMTg3NSBsIDk4LjczMzc2MDEsNC44NDkyMTg3NSBsIDk5LjIxODM4NzUsNC44NDkyMTg3NSBsIDk5LjcwMTk3NTgsNC44NDkyMTg3NSBsIDEwMC4xNzg1Myw0Ljg0OTIxODc1IGwgMTAwLjYwNDA4OCw0Ljg0OTIxODc1IGwgMTAxLjQwMzgwNyw0Ljg0OTIxODc1IGwgMTAyLjAwMjY1Nyw0Ljg0OTIxODc1IGwgMTAyLjM5NTU2NCw0Ljg0OTIxODc1IGwgMTAzLjA4NjM4NSw1LjA3OTQ5MjE5IGwgMTAzLjQ3MDE3NCw1LjU0MDAzOTA2IGwgMTAzLjkzMDcyMSw1LjU0MDAzOTA2IGwgMTA0LjMxNDUxLDUuMDc5NDkyMTkgbCAxMDQuNzc1MDU3LDQuNjk1NzAzMTMgbCAxMDUuMDA1MzMsNC4wODE2NDA2MyBsIDEwNS4xNTg4NDYsMy4zOTA4MjAzMSBsIDEwNS42MTkzOTIsMy4wMDcwMzEyNSBsIDEwNS44NDk2NjYsMi4zOTI5Njg3NSBsIDEwNi4wMDMxODIsMS43MDIxNDg0NCBsIDEwNi42MTcyNDQsMS40NzE4NzUgbCAxMDcuMzA4MDY0LDEuMzE4MzU5MzggbCAxMDcuNTM4MzM4LDAuNzA0Mjk2ODc1IGwgMTA3LjY5MTg1MywwLjAxMzQ3NjU2MjUgbCAxMDguMTUyNCwwLjAxMzQ3NjU2MjUgbCAxMDguMzgyNjc0LDAuNzA0Mjk2ODc1IGwgMTA4LjUzNjE4OSwxLjMxODM1OTM4IGwgMTA4Ljk5NjczNiwxLjcwMjE0ODQ0IGwgMTA5LjM4MDUyNSwyLjE2MjY5NTMxIGwgMTA5Ljg0MTA3MiwyLjU0NjQ4NDM4IGwgMTEwLjA3MTM0NiwzLjIzNzMwNDY5IGwgMTEwLjA3MTM0Niw0LjA4MTY0MDYzIGwgMTEwLjIyNDg2MSw0LjY5NTcwMzEzIGwgMTEwLjY4NTQwOCw1LjA3OTQ5MjE5IGwgMTExLjA2OTE5Nyw1LjU0MDAzOTA2IGwgMTExLjY4MzI2LDUuNjkzNTU0NjkgbCAxMTIuMzc0MDgsNS45MjM4MjgxMyBsIDExMi43NTc4NjksNi4zODQzNzUgbCAxMTMuMjE4NDE2LDYuMzg0Mzc1IGwgMTEzLjYwMjIwNSw1LjkyMzgyODEzIGwgMTE0LjIxNjI2Nyw1LjY5MzU1NDY5IGwgMTE1LjAwNDYyNyw1LjY5MzU1NDY5IGwgMTE1LjkwNDkzOSw1LjY5MzU1NDY5IGwgMTE2LjU5NTc2LDUuOTIzODI4MTMgbCAxMTYuOTc5NTQ5LDYuMzg0Mzc1IGwgMTE3LjU5MzYxMSw2LjUzNzg5MDYzIGwgMTE4LjI1MTc4Miw2LjUzNzg5MDYzIGwgMTE4LjczMjA5Myw2LjUzNzg5MDYzIGwgMTE5LjM4OTA1NCw2LjUzNzg5MDYzIGMgMTE5LjYxMDA5Miw2LjUzNzg5MDYzIDExOS44MzExMyw2LjUzNzg5MDYzIDEyMC4wNTIxNjksNi41Mzc4OTA2MyBjIDEyMC4yMDM3NjUsNi41Mzc4OTA2MyAxMjAuMzgxODAzLDYuNTM3ODkwNjMgMTIwLjU3MjQ0Niw2LjUzNzg5MDYzIGMgMTIwLjc2NDkyMiw2LjUzNzg5MDYzIDEyMC45NTczOTgsNi41Mzc4OTA2MyAxMjEuMTQ5ODc0LDYuNTM3ODkwNjMgYyAxMjEuMzUyMTM1LDYuNTM3ODkwNjMgMTIxLjYzNTY4Niw2LjUzNzg5MDYzIDEyMS45MDI1MTMsNi41Mzc4OTA2MyBsIDEyMi42NTk2MjcsNi41Mzc4OTA2MyBjIDEyMi44OTc1Nyw2LjUzNzg5MDYzIDEyMy4zNTA0NDcsNi43NjgxNjQwNiAxMjMuMzUwNDQ3LDYuNzY4MTY0MDYgbCAxMjMuNzM0MjM2LDcuMjI4NzEwOTQgbCAxMjQuMTk0NzgzLDcuNjEyNSBsIDEyNC41Nzg1NzIsOC4wNzMwNDY4OCBsIDEyNS4wMzkxMTksOC40NTY4MzU5NCBsIDEyNS4yNjkzOTIsOS4xNDc2NTYyNSBsIDEyNS4wMzkxMTksOS43NjE3MTg3NSBsIDEyNC41Nzg1NzIsMTAuMTQ1NTA3OCBsIDEyNC4xOTQ3ODMsMTAuNjA2MDU0NyBsIDEyMy43MzQyMzYsMTAuOTg5ODQzOCBsIDEyMy4zNTA0NDcsMTEuNDUwMzkwNiBsIDEyMi42NTk2MjcsMTEuNjAzOTA2MyBsIDEyMS45ODgyLDExLjYwMzkwNjMgbCAxMjEuMzQ5Mzg0LDExLjYwMzkwNjMgbCAxMjAuNzgwOTg4LDExLjYwMzkwNjMgYyAxMjAuNTkzNjAxLDExLjYwMzkwNjMgMTIwLjQwNjIxNCwxMS42MDM5MDYzIDEyMC4yMTg4MjcsMTEuNjAzOTA2MyBsIDExOS40NDI0NDgsMTEuNjAzOTA2MyBsIDExOC43NzQxMzcsMTEuNjAzOTA2MyBsIDExOC4xNzYyNDYsMTEuNjAzOTA2MyBsIDExNy41MjAyNDUsMTEuNjAzOTA2MyBsIDExNi44MzY2NjcsMTEuNjAzOTA2MyBsIDExNi4yMzg3NzYsMTEuNjAzOTA2MyBjIDExNi4wNTE4NjMsMTEuNjAzOTA2MyAxMTUuODY4NjQ5LDExLjYwMzkwNjMgMTE1LjY5MDkyMiwxMS42MDM5MDYzIGMgMTE1LjQ5MDYyMiwxMS42MDM5MDYzIDExNS4yOTcyOTIsMTEuNjAzOTA2MyAxMTUuMTEzNDk0LDExLjYwMzkwNjMgbCAxMTQuNTM2MDY2LDExLjYwMzkwNjMgYyAxMTQuMzc0OTA3LDExLjYwMzkwNjMgMTE0LjIyNjQ0MSwxMS42MDM5MDYzIDExNC4wOTMzMjIsMTEuNjAzOTA2MyBjIDExMy42NDU5NjcsMTEuNjAzOTA2MyAxMTMuMzcxOTMyLDExLjYwMzkwNjMgMTEzLjM3MTkzMiwxMS42MDM5MDYzIGwgMTEyLjc1Nzg2OSwxMS44MzQxNzk3IGwgMTEyLjM3NDA4LDEyLjI5NDcyNjYgbCAxMTEuOTEzNTMzLDEyLjY3ODUxNTYgbCAxMTEuNTI5NzQ0LDEzLjEzOTA2MjUgbCAxMTEuMDY5MTk3LDEzLjUyMjg1MTYgbCAxMTAuNjg1NDA4LDEzLjk4MzM5ODQgbCAxMTAuMjI0ODYxLDE0LjM2NzE4NzUgbCAxMTAuMDcxMzQ2LDE1LjA1ODAwNzggbCAxMTAuMDcxMzQ2LDE1Ljk5MDI4OTYgbCAxMTAuMDcxMzQ2LDE2Ljc0NjY3OTcgYyAxMTAuMDcxMzQ1LDE3LjIxODkyMjggMTEwLjIyNDg2MSwxNy4zNjA3NDIyIDExMC4yMjQ4NjEsMTcuMzYwNzQyMiBsIDExMC42ODU0MDgsMTcuNzQ0NTMxMyBsIDExMC42ODU0MDgsMTguMjA1MDc4MSBsIDExMC4yMjQ4NjEsMTguNTg4ODY3MiBsIDExMC4wNzEzNDYsMTkuMjc5Njg3NSBsIDExMC4wNzEzNDYsMjAuMTczMDQ4MyBjIDExMC4wNzEzNDYsMjAuMzk1MDU1MiAxMTAuMDcxMzQ2LDIwLjYzNDY3NjYgMTEwLjA3MTM0NiwyMC44ODc3MTk4IGwgMTEwLjA3MTM0NiwyMS42NTY5ODQ4IGwgMTEwLjA3MTM0NiwyMi4zOTA1MjAyIGwgMTEwLjA3MTM0NiwyMy4xODczNjE4IGwgMTEwLjA3MTM0NiwyMy45NzA5MzQ2IGwgMTEwLjA3MTM0NiwyNC42NTM1NTM0IGwgMTEwLjA3MTM0NiwyNS40NzM3MzUxIGwgMTEwLjA3MTM0NiwyNi4yNjUzODEgbCAxMTAuMDcxMzQ2LDI2LjkxMzIyOTQgbCAxMTAuMDcxMzQ2LDI3LjU0NDI5MjEgYyAxMTAuMDcxMzQ2LDI4LjE3NDI1NDIgMTEwLjA3MTM0NiwyOC41NjczODI4IDExMC4wNzEzNDYsMjguNTY3MzgyOCBsIDExMC4yMjQ4NjEsMjkuMTgxNDQ1MyBsIDExMC42ODU0MDgsMjkuNTY1MjM0NCBsIDExMC45MTU2ODIsMzAuMjU2MDU0NyBsIDExMC45MTU2ODIsMzEuMTAwMzkwNiBsIDExMS4wNjkxOTcsMzEuNzE0NDUzMSBsIDExMS41Mjk3NDQsMzIuMDk4MjQyMiBsIDExMS43NjAwMTcsMzIuNzg5MDYyNSBsIDExMS43NjAwMTcsMzMuNjMzMzk4NCBsIDExMS45MTM1MzMsMzQuMjQ3NDYwOSBsIDExMi4zNzQwOCwzNC42MzEyNSBsIDExMi42MDQzNTMsMzUuMzIyMDcwMyBsIDExMi43NTc4NjksMzUuOTM2MTMyOCBsIDExMy4zNzE5MzIsMzYuMDg5NjQ4NCBsIDExNC4wNjI3NTIsMzYuMzE5OTIxOSBsIDExNC40NDY1NDEsMzYuNzgwNDY4OCBsIDExNS4wNjA2MDMsMzYuOTMzOTg0NCBsIDExNS43NDIzMTgsMzYuOTMzOTg0NCBsIDExNi4yMjY4NjYsMzYuOTMzOTg0NCBsIDExNi42ODM5MTcsMzYuOTMzOTg0NCBsIDExNy4wNDE4NTIsMzYuOTMzOTg0NCBsIDExNy43MjczNDgsMzYuOTMzOTg0NCBsIDExOC40Mzc5NDcsMzYuOTMzOTg0NCBsIDExOS4xMjg3NjcsMzYuNzgwNDY4OCBsIDExOS41MTI1NTcsMzYuMzE5OTIxOSBsIDExOS45NzMxMDMsMzUuOTM2MTMyOCBsIDEyMC4zNTY4OTIsMzUuNDc1NTg1OSBsIDEyMC44MTc0MzksMzUuMDkxNzk2OSBsIDEyMS4wNDc3MTMsMzQuNDc3NzM0NCBsIDEyMS4yMDEyMjgsMzMuNzg2OTE0MSBsIDEyMS42NjE3NzUsMzMuNDAzMTI1IGwgMTIxLjg5MjA0OSwzMi43ODkwNjI1IGwgMTIxLjg5MjA0OSwzMS45NDQ3MjY2IGwgMTIyLjA0NTU2NCwzMS4yNTM5MDYzIGwgMTIyLjUwNjExMSwzMC44NzAxMTcyIGwgMTIyLjczNjM4NSwzMC4yNTYwNTQ3IGwgMTIyLjg4OTksMjkuNTY1MjM0NCBsIDEyMy4zNTA0NDcsMjkuMTgxNDQ1MyBsIDEyMy43MzQyMzYsMjguNzIwODk4NCBsIDEyNC4zNDgyOTksMjguNDkwNjI1IGwgMTI1LjAzOTExOSwyOC4zMzcxMDk0IGwgMTI1LjQyMjkwOCwyNy44NzY1NjI1IGwgMTI2LjAzNjk3MSwyNy42NDYyODkxIGwgMTI2Ljg4MTMwNywyNy42NDYyODkxIGwgMTI3LjU3MjEyNywyNy44NzY1NjI1IGwgMTI3Ljk1NTkxNiwyOC4zMzcxMDk0IGwgMTI4LjU2OTk3OCwyOC40OTA2MjUgbCAxMjkuMjYwNzk5LDI4LjcyMDg5ODQgbCAxMjkuNjQ0NTg4LDI5LjE4MTQ0NTMgbCAxMzAuMTA1MTM1LDI5LjU2NTIzNDQgbCAxMzAuMzM1NDA4LDMwLjI1NjA1NDcgbCAxMzAuNDg4OTI0LDMwLjg3MDExNzIgbCAxMzAuOTQ5NDcxLDMxLjI1MzkwNjMgbCAxMzAuOTQ5NDcxLDMxLjcxNDQ1MzEgbCAxMzAuNDg4OTI0LDMyLjA5ODI0MjIgeiBcIjtcbnZhciBsZXR0cmVTID0gXCJtIDEzOC44MTE1MTcsMzUuMTA2MTAxMSBjIDEzOC42Mzc4MSwzNS43MzY1MTExIDEzOC42MDY1NjQsMzUuODAxMDYxNCAxMzguMTk4NDU1LDM2LjIwOTE3MDQgYyAxMzguMDU2NzUzLDM2LjM1MDg3MjIgMTM3Ljk5OTE1NSwzNy4wMjE1OTI4IDEzNy45NDg3OTgsMzcuMjIzMDE5OSBjIDEzNy44Njk1MzcsMzcuNTQwMDYyMiAxMzcuOTI2NDExLDM3LjgyOTIzMzEgMTM3Ljk0ODc5OCwzOC4wNzc1NDc1IGMgMTM3Ljk4MzA3OCwzOC40NTc3ODA2IDEzNy45Mjc3NjYsMzguOTQxNTM2NyAxMzguMTI2MDM3LDM5LjIyMTc0NzkgYyAxMzguMjQwNzEzLDM5LjM4MzgxNjQgMTM4LjQzODU5NiwzOS40MzM4MzUgMTM4LjU2MDU0NCwzOS42MjAxOTY0IGMgMTM4Ljc0NzA2NiwzOS45MDUyNDA4IDEzOC43OTU5MSw0MC40NTg0NTI3IDEzOC45NTE1OTksNDAuODMzOTMzIGMgMTM5LjA2MTEzMSw0MS4wOTgwOTIyIDEzOS4yNjM5MTgsNDEuMjA0NTk1OCAxMzkuNDI5NTU2LDQxLjM3MDIzNDMgYyAxMzkuNTkyNjcyLDQxLjUzMzM0OTUgMTM5LjYyNjQxMiw0MS45MzQ1MzYzIDEzOS42MjIxNjksNDIuNDI0NTg3NiBjIDEzOS42MTgzMiw0Mi44NjkwNjI1IDEzOS41ODMyMjYsNDMuMzg2NjQxIDEzOS41ODUwNzYsNDMuODY1OTk0NyBjIDEzOS41ODY1NzksNDQuMjU1NDU5NCAxMzkuNjEyNDcxLDQ0LjYxOTY5MTIgMTM5LjY5OTMyNCw0NC44OTg5ODA2IGMgMTM5Ljc5ODUxOCw0NS4yMTc5NTU3IDE0MC4wNTg1MzMsNDUuMzM2ODUxNSAxNDAuMjQ5MzgyLDQ1LjUxNjI0MzggYyAxNDAuNDY4MzcyLDQ1LjcyMjA4NzcgMTQwLjY2NzY5OCw0NS45NTI3NDA1IDE0MC44ODM3MTQsNDYuMjA2NDE5NiBjIDE0MS4yODc5ODcsNDYuNjgxMTc5MSAxNDEuMzYzNTYzLDQ2LjgwNDI2NjIgMTQxLjc1ODczMyw0Ni44ODg5NDUyIGMgMTQyLjE1MzkwMyw0Ni45NzM2MjQxIDE0Mi4zODk1NzEsNDYuODg4OTQ1MiAxNDIuNTc3Myw0Ni44ODg5NDUyIGMgMTQyLjkxOTAxOCw0Ni44ODg5NDUyIDE0My4wNjI1ODYsNDYuNTQ0NzkyNCAxNDMuMTcwMDU1LDQ2LjQzNzMyMjkgYyAxNDMuNTI0NCw0Ni4wODI5Nzc2IDE0NC4zMDI2MSw0Ni4wODQ5MzI3IDE0NC40NTYxMDcsNDUuOTU3NjkyOCBjIDE0NC42MDk2MDMsNDUuODMwNDUyOCAxNDQuNTM4MjUyLDQ1Ljc4ODg5OCAxNDQuNjgwMTY5LDQ1LjY0Njk4MTEgYyAxNDQuODY3MzkxLDQ1LjQ1OTc1OTIgMTQ1LjM0MDgsNDUuMzQ5MjM5NSAxNDUuOTIxMDQ3LDQ1LjI5MTkzMjkgYyAxNDYuNDE2ODc5LDQ1LjI0Mjk2MzMgMTQ2Ljk5MDcyMyw0NS4yMzI4NSAxNDcuNTMwNjczLDQ1LjI0NjkzNjUgYyAxNDguMjM0NTk5LDQ1LjI2NTMwMDggMTQ4Ljg4MDkxNyw0NS4zMjQ3OTUyIDE0OS4yMjE2NjEsNDUuMzkyOTQ0MiBjIDE0OS43OTM1MjcsNDUuNTA3MzE3MiAxNDkuNTkyMjkzLDQ1Ljk0MDY3OTQgMTUwLjMxMTM0NCw0Ni4wODQ0ODk2IGMgMTUwLjU5MDYxNSw0Ni4xNDAzNDM3IDE1MC45ODM1MDgsNDYuMTMyODU5NiAxNTEuMzg1NzA1LDQ2LjEzMzk4MzMgYyAxNTEuNzQwNDYxLDQ2LjEzNDk3NDUgMTUyLjEwMjQ1NSw0Ni4xNDI2NjI2IDE1Mi40MDAxMDEsNDYuMjA2NDE5NiBjIDE1My4wMzUxOTUsNDYuMzQyNDU5NiAxNTIuODQzMjU5LDQ2LjY5ODk2NDkgMTUzLjMwMzc5Nyw0Ni44NTUzNTYzIGMgMTUzLjc0MDEwOSw0Ny4wMDM1MjA4IDE1NC4xMjY2NTksNDcuMDYwNzE3NyAxNTQuNDYyNzY3LDQ3LjA1MjUzMDcgYyAxNTQuODQ4ODU4LDQ3LjA0MzEyNjMgMTU1LjE2ODM4OSw0Ni45NDc0NDUxIDE1NS40MjAzMjksNDYuODA0MjY2MiBjIDE1NS43MDExNTMsNDYuNjQ0NjcyMyAxNTUuODU2NTgyLDQ2LjMwMTg4NTIgMTU2LjIzODQ0NCw0Ni4yMDY0MTk2IGMgMTU2LjYzNzU1LDQ2LjEwNjY0MyAxNTcuMDYyNDYzLDQ2LjA3NzM4NjQgMTU3LjUwMTE5Myw0Ni4wNzk5MDAzIGMgMTU3Ljk2MjUxNCw0Ni4wODI1NDM2IDE1OC40MzkxMTEsNDYuMTIwMzEzMyAxNTguOTE3MDQ3LDQ2LjE0ODE2MDkgYyAxNTkuMzk1MTU4LDQ2LjE3NjAxODcgMTU5Ljg3NDYwOSw0Ni4xOTM5NDcgMTYwLjM0MTQ0NSw0Ni4xNTY4NDc3IGMgMTYwLjY5NTYyNCw0Ni4xMjg3MDExIDE2MS4wNDI1NDMsNDYuMDY4ODgwOSAxNjEuMzc2MTA3LDQ1Ljk1NzY5MjggYyAxNjEuNTU2MTY1LDQ1Ljg5NzY3MzQgMTYxLjU5NTUwMSw0NS42NTQyMTU5IDE2MS44NTU5NTcsNDUuNTE2MjQzOCBjIDE2Mi4yMTQzNTcsNDUuMzI2Mzg3MiAxNjIuNzYwNDUsNDUuMjA0NTcyMiAxNjIuOTk5MTI3LDQ1LjEyNTAxMzMgYyAxNjMuMzAxNDY5LDQ1LjAyNDIzMjYgMTYzLjI1NDIxMSw0NC44MTc2ODQzIDE2My41NDI3OTgsNDQuNjczMzkxIGMgMTYzLjgxMTkyMSw0NC41Mzg4MjkzIDE2NC4xOTg3NTEsNDQuNTA0MDI5NCAxNjQuNTM3NDY5LDQ0LjM2Mjg5OTEgYyAxNjQuODc2MTg3LDQ0LjIyMTc2ODcgMTY1LjI0NTYyOCw0My43OTU4NzE4IDE2NS41OTU5Niw0My40NDU1NDA3IGMgMTY1LjY5NDA3NCw0My4zNDc0MjYzIDE2Ni4wODk0MjcsNDIuOTUxNTc3IDE2Ni40MTQwMzEsNDIuNTg0NjMzNyBjIDE2Ni43Mzg2MzUsNDIuMjE3NjkwMyAxNjYuNjQ5Mzk4LDQxLjY5NTgzODEgMTY2LjgyMzgxLDQxLjM3MDIzNDMgYyAxNjYuOTQyNTg3LDQxLjE0ODQ5MjggMTY3LjIzMTcyOCw0MS4wNzE3MjAyIDE2Ny4zNDYxNjQsNDAuODMzOTMzIGMgMTY3LjUxNDc5Nyw0MC40ODM1Mjg2IDE2Ny40NzA1NjQsMzkuOTU4NDMzOSAxNjcuNjQyNTQyLDM5LjYyMDE5NjQgYyAxNjcuNzUyMzMyLDM5LjQwNDI2ODEgMTY4LjA2NDc1NiwzOS4zMzA0NjQ2IDE2OC4xNjQ3MzIsMzkuMTI2MjMyOCBjIDE2OC4zMzc1NzYsMzguNzczMTQzMiAxNjguMzgyMjM1LDM4LjIxOTI0MjUgMTY4LjUxNzM5NSwzNy45NTQ4Mzc1IGMgMTY4LjY2MzgwMiwzNy42Njg0Mjk0IDE2OC45NDYzNTMsMzcuNTM2Njk1MSAxNjkuMDExMzU4LDM3LjQwNjY4NTcgYyAxNjkuMTYxOTksMzcuMTA1NDIxNCAxNjkuMTk2NDE4LDM2LjU4MzYxOTQgMTY5LjIxNjAwMSwzNi4wNTM0MzAxIGMgMTY5LjIzNjgzLDM1LjQ4OTUwNTUgMTY5LjI0MDg2NiwzNC45MTYwOTI0IDE2OS4zNTAwNzQsMzQuNTg4NDY5MSBjIDE2OS40Njk3NzQsMzQuMjI5MzY5MyAxNjkuODE1ODEsMzQuMjkwNjkwMiAxNjkuOTcxMDU4LDMzLjc2ODUwMjcgYyAxNzAuMTI2MzA2LDMzLjI0NjMxNTIgMTY5LjY3ODYxMywzMi43Njc2MzU0IDE2OS40MDkyOSwzMi40OTgzMTI5IGMgMTY5LjIwOTMwNSwzMi4yOTgzMjczIDE2OS4zNDEwMzUsMzEuNDQzNzY1MyAxNjguNTczODUsMzAuODQ3MDY1OSBjIDE2OC4zMTk4MTEsMzAuNjQ5NDgwNSAxNjguNDU2NTU1LDMwLjAwODY2NjcgMTY4LjE2NDczMiwyOS42MDUzMjM0IGMgMTY3Ljk5MzQzNywyOS4zNjg1NjkxIDE2Ny43MjQxNjksMjkuMjg1NjczOSAxNjcuNjQyNTQyLDI5LjA0MDc5NDYgYyAxNjcuNDU1ODk2LDI4LjQ4MDg1NjggMTY3LjYwOTI3NSwyOC40Mjc5NjIyIDE2Ny4zNDYxNjQsMjcuOTgyMzA0IGMgMTY3LjEyNDE5NCwyNy42MDYzMjkzIDE2Ni42MjI3OCwyNy4xNDU5NzU5IDE2Ni4yODc2NzMsMjYuODU1NDA4OCBjIDE2NS45NDQ2MzUsMjYuNTU3OTY0MyAxNjUuNzExMDU0LDI2LjI0ODQ5NDcgMTY1LjQ1NDMzMSwyNi4wMDY0NTMxIGMgMTY1LjE0NDYyOCwyNS43MTQ0NjEzIDE2NS4wNDUwNjEsMjUuNDg0NTk0IDE2NC42Nzg3NjYsMjUuMjcxMTAzNCBjIDE2NC4zNDIyNywyNS4wNzQ5ODA1IDE2My45MjIzOTgsMjUuMDQ2NjYxNyAxNjMuNTQyNzk4LDI0Ljg3NzM5NTUgYyAxNjMuMTk0NDc4LDI0LjcyMjA3NzMgMTYzLjA4Mzk3MiwyNC4zNTUyMDYyIDE2Mi41MDc1MzYsMjQuMjI4MTg3NyBjIDE2MS45MzExLDI0LjEwMTE2OTMgMTYxLjQ4NjUwNSwyNC4yMjc4OTM2IDE2MS4wMjM0NDIsMjQuMDcyOTQxOCBjIDE2MC42ODM5MjYsMjMuOTU5MzMxMyAxNjAuNjEzNzIyLDIzLjU2MDAxNjkgMTYwLjI0NzIxNCwyMy40Mzc4NDc4IGMgMTU5Ljc4ODAzMiwyMy4yODQ3ODcgMTU5LjE0Mjg3OSwyMy4zMDc5MDMzIDE1OC40NzAxMjIsMjMuMzQ3OTExNSBjIDE1Ny45NDQ2MTIsMjMuMzc5MTYzIDE1Ny40MDIyNiwyMy40MjA3MjEzIDE1Ni45MTg1NDUsMjMuMzk2NjY4MyBjIDE1Ni41NDA1MDcsMjMuMzc3ODcwMSAxNTYuMTk4Mjg0LDIzLjMxODk5NyAxNTUuOTI3OTA4LDIzLjE4MzgwOTEgYyAxNTUuNzA2NjU2LDIzLjA3MzE4MjkgMTU1LjQ5MTA2MSwyMi42MjU0MzA4IDE1NC45MTI0MTgsMjIuNTQ4NzEzMyBjIDE1NC4zMzM3NzYsMjIuNDcxOTk1NyAxNTMuODE1MjY5LDIyLjU0MDg4OTkgMTUzLjQ3Nzc0MywyMi4zOTUzNzM0IGMgMTUzLjE2Mjg1NywyMi4yNTk2MTc3IDE1My4wNTk3NjMsMjEuODY2OTYxMiAxNTIuOTA3NjgsMjEuODI4OTQwMyBjIDE1Mi42MTk4ODUsMjEuNzU2OTkxNSAxNTIuMzQ0Mjk2LDIxLjcwODk4MzcgMTUyLjA3OTcxMSwyMS42Nzc3NjcxIGMgMTUxLjcyMzQ4OSwyMS42MzU3Mzg2IDE1MS4zODcyMTMsMjEuNjI0MTQ2NyAxNTEuMDY3OTUxLDIxLjYyNTU0MjIgYyAxNTAuNTM5NzEyLDIxLjYyNzg1MTIgMTUwLjA1ODA1LDIxLjY2NTcxNDMgMTQ5LjYwOTY3NSwyMS42NjAwOTYxIGMgMTQ5LjE2NDI1NiwyMS42NTQ1MTQ5IDE0OC43NTE2ODgsMjEuNjA2MDI0IDE0OC4zNTg5NDIsMjEuNDM3MTQxMyBjIDE0Ny45NzU5ODcsMjEuMjcyNDY5IDE0Ny44OTg4MjYsMjEuMDk1MDUyNyAxNDcuMjA4NjUsMjAuNDE3NjE4MyBjIDE0Ni45OTY0NzYsMjAuMjA5MzYxNCAxNDYuNzk5MTc1LDIwLjAyNjQ0NjYgMTQ2LjYyMDQ4MiwxOS44NjIzMTMzIGMgMTQ2LjIxNzkwOSwxOS40OTI1NDEzIDE0NS45MDk3ODQsMTkuMjE4MDk0NCAxNDUuNzM4ODI2LDE4Ljk2Mzk1NjggYyAxNDUuNDkxOTg0LDE4LjU5NzAxMzQgMTQ1LjU4MzU4MSwxOC4xNDUzOTAyIDE0NS4zNjgxODIsMTcuNzkyNTU5NyBjIDE0NS4xNTI3ODMsMTcuNDM5NzI5MiAxNDUuMTY1OTI2LDE3Ljc5ODQ5NTMgMTQ0Ljg5NjA5NCwxNi45ODg5OTk1IGMgMTQ0LjgxMDU5MywxNi43MzI0OTUyIDE0NS4yODc1OTksMTYuNjAxNDQzNCAxNDUuMzY4MTgyLDE2LjI5ODgyMzcgYyAxNDUuNjM5MzcyLDE1LjI4MDQwNjcgMTQ1LjY0Mzk1NywxNS4xODUwOTY2IDE0Ni4xNjU5ODIsMTQuNjU5NDI0NiBjIDE0Ni4zNjAxMjcsMTQuNDYzOTIyOCAxNDYuNDkwNjI5LDEzLjQ5NzY3NjUgMTQ2LjcyNjA5LDEzLjI2MjIxNjQgYyAxNDYuOTU4NzM5LDEzLjAyOTU2NzQgMTQ3LjIyMjQzOSwxMi44MjMwOTIxIDE0Ny40NTkyNjksMTIuNTkwNTU0MyBjIDE0Ny42OTE5MDgsMTIuMzYyMTMxNCAxNDcuODk4NjIxLDEyLjEwODU2MDEgMTQ4LjAyNDUwNiwxMS43ODAzMjgzIGMgMTQ4LjI3ODU0NSwxMS4xMTc5NTE1IDE0Ny45ODIxNjcsMTEuMDc5MTA1MyAxNDguMzU4OTQyLDEwLjkwNTMwODkgYyAxNDguNzM1NzE3LDEwLjczMTUxMjUgMTQ4Ljc3MTMwOSwxMS4xNjA5NDIgMTQ5LjIyMTY2MSwxMS4zOTkyNzE2IGMgMTQ5LjUwNDU2NywxMS41NDg5ODczIDE0OS44MjIzMTQsMTEuNTIwMzY4NyAxNTAuMTk3OTQxLDExLjU1NDUxNjYgYyAxNTAuNTQ2NDEsMTEuNTg2MTk1NiAxNTEuMTI5Mjk3LDExLjU3Nzg3ODcgMTUxLjcxNjM5NywxMS41Njc2MjgxIGMgMTUyLjYwNjAyNCwxMS41NTIwOTU2IDE1My41MDUzMjIsMTEuNTMyMTIzNCAxNTMuNjEzMzQsMTEuNjQwMTQwOSBjIDE1My43NjY1MDgsMTEuNzkzMzA4NSAxNTMuOTcxNTc3LDEyLjA4NjQ2NTQgMTU0LjIwNjA5NiwxMi4yMDM3MjQ5IGMgMTU0LjY5ODMxNSwxMi40NDk4MzQ0IDE1NS4zMDY3MywxMi40MjY2NjA4IDE1NS45MjAwNzMsMTIuMzg1NTkzMyBjIDE1Ni42MTQ3NTUsMTIuMzM5MDc5NyAxNTcuMzE1NzU5LDEyLjI2OTYxMTUgMTU3Ljg2MTQyLDEyLjU0MjQ0MjEgYyAxNTcuOTg0NzI4LDEyLjYwNDA5NjIgMTU4LjIzMzM3MywxMi45NTgzNDkzIDE1OC4zNDE0NjgsMTMuMDEyMzk2OSBjIDE1OS4wMjQyODQsMTMuMzUzODA0OCAxNTkuMjc4MjA0LDEzLjI2MTU3ODggMTU5LjYzOTY4NSwxMy40NTg4NzU4IGMgMTU5Ljk0MjMwOCwxMy42MjQwNDczIDE2MC4xMTM1NTksMTMuODk2NTM2NCAxNjAuNDg2NDc4LDE0LjI3MTA5MjcgYyAxNjAuNjg2ODAxLDE0LjQ3MjI5NDYgMTYwLjg4NjY1MiwxNC42ODU5ODEgMTYxLjE4NDcsMTQuOTg0MDI4NSBjIDE2MS4zNjU3ODMsMTUuMTY1MTExMyAxNjEuNjk1NzY1LDE1LjUzNzA3NjggMTYyLjA1ODk1NywxNS45MzI2Mjg3IGMgMTYyLjI2NzM5LDE2LjE1OTYzMzYgMTYyLjQ4Njc2MSwxNi4zOTQ0MDY3IDE2Mi42OTUyMDIsMTYuNjA1MzI2NyBjIDE2Mi45ODE1MywxNi44OTUwNjA4IDE2My4yNDcyMzQsMTcuMTM5Nzg0OCAxNjMuNDM1NjMzLDE3LjI1NzUzMzggYyAxNjQuMDAwMTYyLDE3LjYxMDM2NDMgMTY0LjAyMTI4MiwxNy4zMjUyMDUgMTY0LjMzOTM4NywxNy40MTE1MDI2IGMgMTY0LjgxMDIzOCwxNy41MzkyMzc4IDE2NS4yNjQwNDcsMTcuNDQxNDYyMyAxNjUuNjU3MTIzLDE3LjE5ODQ4MDUgYyAxNjUuOTk1NjcsMTYuOTg5MjA2MyAxNjYuMjg5MTY2LDE2LjY3MjIxODkgMTY2LjUwOTY5NywxNi4yOTg4MjM3IGMgMTY2LjcwMTk1OCwxNS45NzMyOTUgMTY2Ljc0MTAzLDE1LjQ4ODA1NTEgMTY2LjczNDEwNCwxNC45MTM5OTI3IGMgMTY2LjcyODA4NCwxNC40MTQ5OTM4IDE2Ni42ODczMSwxMy44NDg4ODIyIDE2Ni42ODIxODIsMTMuMjYyMjE2NCBjIDE2Ni42NzkwMzYsMTIuOTAyMzg0OCAxNjYuNjk1MTk3LDEyLjUwMDY4OTggMTY2LjY4NzY5NSwxMi4xMjk0MzI1IGMgMTY2LjY3NzQzMiwxMS42MjE1NDA1IDE2Ni42MjI4ODIsMTEuMTcwNjEzIDE2Ni40MTQwMzEsMTAuOTYxNzYxNyBjIDE2Ni4xMjgwMjgsMTAuNjc1NzU5MiAxNjYuMTU5MDc2LDEwLjc5MDU3MDkgMTY1Ljk2MjQwNywxMC4zOTcyMzI5IGMgMTY1Ljg1NjIxMiwxMC4xODQ4NDQ0IDE2NS45MjQ4NzcsOS41OTI0MjIwNiAxNjUuNDU0MzMxLDkuMTY5MzgyOTcgYyAxNjQuOTgzNzg1LDguNzQ2MzQzODggMTY0Ljk3NDQ4MSw4LjM2MTgwNzk3IDE2NC40NTIyOTIsOC4yNjYxMzY4NSBjIDE2NC4wNjYzMDcsOC4xOTU0MjAwMiAxNjMuNjg1MjIzLDguMTQ5MDk4NTIgMTYzLjMwNjM5Nyw4LjEyMDExMjA5IGMgMTYyLjc5MTk3Myw4LjA4MDc1MDIzIDE2Mi4yODE3MTEsOC4wNzMzNTQzNSAxNjEuNzY4OTk0LDguMDgwMjQ1MjggYyAxNjEuMjk2MDI2LDguMDg2NjAxOTYgMTYwLjgyMDk2OSw4LjEwNTExNjA3IDE2MC4zMzg2MjYsOC4xMjE5MDk3NyBjIDE1OS44NTE0MTcsOC4xMzg4NzI5MSAxNTkuMzU2Nzc1LDguMTU0MDgwNzYgMTU4Ljg0OTM0NSw4LjE1MzIzMTE2IGMgMTU4LjMwNDcxLDguMTUyMzE5MjcgMTU3LjUyMzUsOC4xOTc5MjQ3NSAxNTcuMDI4NzQsOC4wNzQyMzQ4OCBjIDE1Ni40Nzg2NzMsNy45MzY3MTc5NCAxNTYuNDEwNzU3LDcuNTM3OTUzNTIgMTU1LjkyNzkwOCw3LjM3NzAwNDEgYyAxNTUuNDQ4NjI0LDcuMjE3MjQyNTUgMTU0LjY1OTkwOSw3LjM3NzAwNDEgMTUzLjcwNzgsNy4yNzg4NzMxMyBjIDE1My4yNDExNjksNy4yMzA3Nzg4NSAxNTMuMTc3NDI0LDYuOTEzNjUxOTQgMTUyLjkwNzY4LDYuNzI4NDU3NiBjIDE1Mi42Mzc2ODgsNi41NDMwOTMxNCAxNTIuMjMwNTUsNi40MTc5NDA4NiAxNTEuNTY2OTI0LDYuNDc0NDE5NjYgYyAxNTEuMTIwNzgsNi41MTIzODkzMyAxNTAuODMyMDcsNi41NzMyMTIxOCAxNTAuMzcxOTUzLDcuMTIzNjI3NzEgYyAxNDkuODAyNzcxLDcuMzYzNTUyNDEgMTQ5LjI3OTY5Miw3LjM3MTA1NDMgMTQ5LjAxMjQzMSw3LjUwNDY4NDYxIGMgMTQ4Ljg3OTA4OSw3Ljk0OTE2MTU4IDE0OC40MzQ3NDYsOC4wOTMzNDgxNiAxNDcuOTUzMDE4LDguMTUwNzAxOTMgYyAxNDcuMzY0MDIzLDguMjIwODI2NjMgMTQ2LjcxOTEzOSw4LjE2MTE0MzIyIDE0Ni41MTg0NzUsOC4zNjE4MDc5NiBjIDE0Ni4zOTQ4MjEsOC40ODU0NjE2OSAxNDYuMzE0OTU2LDguNjUzODIyMjYgMTQ2LjE2NTk4Miw4LjgwMjc5NjU5IGMgMTQ1Ljk0ODY3OCw5LjAyMDEwMDg3IDE0NS40NDAxOTEsOC44OTA3OTE1OCAxNDUuMTMwMzA4LDkuMDQ1NzMyNzYgYyAxNDQuODgwNDMsOS4xNzA2NzE5NCAxNDQuNjYwMjYsOS40NDg2Mzc0NiAxNDQuNDU2MTA3LDkuNjUyNzkwNzIgYyAxNDQuMzQ2ODQzLDkuNzYyMDU0MTIgMTQzLjc5ODE3Miw5LjY3MDg0NDkyIDE0My4zMDQzMzUsMTAuMDI5NzMxNiBjIDE0Mi44MTA0OTgsMTAuMzg4NjE4NCAxNDIuMTU3ODY1LDExLjA0MzI1MTIgMTQyLjA4OTg1LDExLjEzNTU1ODIgYyAxNDEuOTU3NjI4LDExLjMxNTAwMTEgMTQxLjY0MzQyNCwxMS41MTA4NzQzIDE0MS41NDM5LDExLjcwOTkyMjQgYyAxNDEuMTAyNzUyLDEyLjU5MjIxNzYgMTQxLjM1NzgyMiwxMi40NzAzODA5IDE0MS4yMzY0OTMsMTIuNjk4MDE0NCBjIDE0MS4xNzY0MDUsMTIuODEwNzUwNSAxNDEuMDYxMzU4LDEzLjA0NzUzOTYgMTQwLjY1MDAyMiwxMy40NTg4NzU5IGMgMTQwLjY0MjQ1NCwxMy40NjY0NDM5IDE0MC40NzUzMjYsMTMuODY3MTQ0NSAxNDAuNDc5NDQxLDE0LjI3MTA5MjQgYyAxNDAuNDgzMjMyLDE0LjY0MzMzOTUgMTQwLjQ5NjI0NCwxNS4wMDU1ODMyIDE0MC40MDk0MjQsMTUuMTc5MjIzNSBjIDE0MC4xMTU0MDEsMTUuNzY3MjY4NiAxNDAuMDc5MjkyLDE1LjUyODg5OTcgMTM5Ljg0ODAwNywxNS45OTM1MzQ0IGMgMTM5Ljc1Nzk0MSwxNi4xNzQ0NzAxIDEzOS41ODAzNiwxNi42ODg5Mjc4IDEzOS42OTkzMjQsMTYuOTg4OTk5MiBjIDEzOS44ODU4NTQsMTcuNDU5NDk3OSAxNDAuMzc0NDc2LDE3LjkwMjg4ODYgMTQwLjQwOTQyNCwxOC4wNzc2Mjc0IGMgMTQwLjQ3NjY4LDE4LjQxMzkwOCAxNDAuNDc5MzM0LDE4LjY2MjgzMyAxNDAuNDc5NDQxLDE5LjExNzAzOTggYyAxNDAuNDc5NTY5LDE5LjY2NjU5NzkgMTQwLjczOTAxOSwxOS45NjIxMzgzIDE0MS4wNTMzNywyMC4yMjY2MjQ1IGMgMTQxLjM0NDY4NiwyMC40NzE3Mjg4IDE0MS42ODMxNTEsMjAuNjkwMTYzNSAxNDEuOTA2MDczLDIxLjA1OTM3OTIgYyAxNDIuMzY5NTQ3LDIxLjgyNzAwNjcgMTQyLjIwNTg4NCwyMi4zMzExOTg2IDE0Mi4zNjk1NDcsMjIuMzk1MzczNCBjIDE0Mi43NjA2MDMsMjIuNTQ4NzEzMyAxNDMuMzI4ODY2LDIyLjQ5NjY1NTYgMTQzLjYwMDY0OSwyMi43Njg0MzggYyAxNDMuODgzODIxLDIzLjA1MTYxIDE0My44MzIzODUsMjMuOTk5NTQwMyAxNDQuMDkzMDg5LDI0LjA3Mjk0MTggYyAxNDQuNzg2NzMsMjQuMjY4MjM2OSAxNDQuOTE4NDUsMjQuMzMyNjYwNCAxNDUuMjY2NDI2LDI0LjQ0ODUyOTUgYyAxNDUuNjE0NDAzLDI0LjU2NDM5ODcgMTQ1LjY4ODAyNSwyNC45MjY0ODc5IDE0Ni4xNjU5ODIsMjUuMzMyMDI2MiBjIDE0Ni42NDM5MzgsMjUuNzM3NTY0NiAxNDYuODg4NTg0LDI2LjAwNjQ1MzEgMTQ3LjIwODY1LDI2LjM3NDg0MTIgYyAxNDcuNDczNTc2LDI2LjY3OTc2MzQgMTQ3LjgyODk2LDI2LjcwMjI2MjYgMTQ4LjE3ODUxMywyNi43OTI2NTA3IGMgMTQ4LjQ1MzI4MywyNi44NjM3MDEgMTQ4LjcyNDQ0OSwyNi45NzY2OTg5IDE0OC45NDUyNDcsMjcuMzAxNzg4OCBjIDE0OS4xMDUyODEsMjcuNTM3NDEzNyAxNDkuNTkyMjQsMjcuNTg5NzYyMiAxNTAuMTc0NzA4LDI3LjU5MTA0MTEgYyAxNTAuNjI1Njg5LDI3LjU5MjAzMTMgMTUxLjEzMzkyNywyNy41NjI0MDYyIDE1MS41OTIwMDcsMjcuNTYzNTMwMiBjIDE1Mi4wNTE2OCwyNy41NjQ2NTgyIDE1Mi40NjA4NDcsMjcuNTk2NzQ5NCAxNTIuNzEwOTcsMjcuNzIxODEwOCBjIDE1My4xMTk3NzYsMjcuOTI2MjEzOCAxNTIuOTQ3ODYsMjcuOTA0NTI0NSAxNTMuMTkwMTY5LDI4LjE0NjgzMzQgYyAxNTMuNTAzNjgyLDI4LjQ2MDM0NiAxNTQuMDczNjU5LDI4LjQ1Njk4MzkgMTU0LjY1NDA2MiwyOC40MzU3MzUzIGMgMTU1LjE4NTQ1NiwyOC40MTYyODEgMTU1LjcyNTU5LDI4LjM4MTgzMzMgMTU2LjA4NTYzNiwyOC41NjE4NTY2IGMgMTU2LjQzNzk2NiwyOC43MzgwMjExIDE1Ni4zNjk4MDUsMjguOTI4NDM2MyAxNTYuOTQwMTY2LDI5LjIxMzYxNjggYyAxNTcuMDc0MzI0LDI5LjI4MDY5NiAxNTcuNDk3MjY0LDI5LjI2NzY5NyAxNTcuOTA4ODQ4LDI5LjI5MjUxODMgYyAxNTguMTg2ODA4LDI5LjMwOTI4MTIgMTU4LjQ1OTU4OCwyOS4zNDMyOTM1IDE1OC42MzQ3NDEsMjkuNDMwODY5NiBjIDE1OS4wNjkyNDgsMjkuNjQ4MTIyNCAxNTkuMDMyNTUsMjkuODU5ODg5NyAxNTkuMTk1NzU2LDI5LjkxNDI5MTkgYyAxNTkuNTc5Nzk3LDMwLjA0MjMwNTQgMTYwLjM4NjYxNiwzMC4xNjY5NjU0IDE2MC42MDk1NTQsMzAuMzg5OTAzMiBjIDE2MC45OTg5NiwzMC43NzkzMDk1IDE2MS4zNTkxODgsMzEuMjUwODYyOCAxNjEuNjQ3MzE5LDMxLjUxNjUwMTMgYyAxNjEuOTg3NTgyLDMxLjgzMDIwMzEgMTYyLjIxMjM4LDMyLjA2MTc1OSAxNjIuMzMyNTM3LDMyLjUxMzUwMyBjIDE2Mi40MDM3MDMsMzIuNzgxMDYyMyAxNjIuNDM4MTYyLDMzLjEyNTg2MjkgMTYyLjQzODE2MiwzMy42MTA3MjA3IGMgMTYyLjQzODE2MiwzNC4wMTg2NjE5IDE2Mi40NTI4MzcsMzQuNzM3NDA4NSAxNjIuMjI2NjYxLDM0Ljk2MzU4NDUgYyAxNjIuMDk5NDE1LDM1LjA5MDgyOTYgMTYxLjM0MzA0LDM1LjE3NDA4NDggMTYxLjE4NDY5OSwzNS4yNTMyNTU1IGMgMTYwLjgwMzE1MywzNS40NDQwMjg3IDE2MC40MzY2MzYsMzUuODk5MDgxNSAxNjAuMTI2NTQ3LDM2LjIwOTE3MDQgYyAxNTkuODQxNzk0LDM2LjQ5MzkyMzMgMTU5LjE3OTE3NywzNy4yNDQyMjUxIDE1OC44NDkzNDQsMzcuNTI3MTcyNyBjIDE1OC41NjIzMjMsMzcuNzczMzk0NyAxNTcuOTQ1ODMsMzcuNjg4NzgzMSAxNTcuNjQ5ODU5LDM3Ljg0NTgxMjkgYyAxNTcuMjg4NTU5LDM4LjAzNzUwMzUgMTU3LjMwNDMzNCwzOC4yMjk0MjMyIDE1Ny4wMjg3MzksMzguMzY3MjIwMyBjIDE1Ni43Nzk2MTUsMzguNDkxNzgyNyAxNTYuMzQ5MjMyLDM4LjUzMjU0NDUgMTU1Ljg2OTkxNSwzOC41NDQ0NjgxIGMgMTU1LjQ4MjQ0MSwzOC41NTQxMDcgMTU1LjA2Mjk4OCwzOC41NDQ5MDA1IDE1NC42ODE0NiwzOC41NDU4ODM4IGMgMTU0LjEzMTg5LDM4LjU0NzMwMDIgMTUzLjY2MTAxLDM4LjU2OTg1OTMgMTUzLjQ3Nzc0MiwzOC43MDAzNDA1IGMgMTUzLjA4OTE3MiwzOC45NzY5OTA0IDE1My4wMTkyMDIsMzkuMjQ2MjkyIDE1Mi40NDg3MjUsMzkuMzM1MDQyNCBjIDE1Mi4xODE1MzksMzkuMzc2NjA5MSAxNTEuODA0NTY0LDM5LjM3ODU3MDggMTUxLjIzMzY0OCwzOS4zMjMxMzM0IGMgMTUwLjQxODc0MSwzOS4yNDQwMDM5IDE1MC41OTYxMDYsMzguODU5NjYwNSAxNTAuMTk3OTQxLDM4LjcwMDM0MDUgYyAxNTAuMDI5MDA1LDM4LjYzMjc0MyAxNDkuNjgwMzg0LDM4LjU5NTM3NDEgMTQ5LjMyMTAzMywzOC41NjA4MTk4IGMgMTQ4LjgzMzQzMywzOC41MTM5MzM0IDE0OC4zMjYwNzksMzguNDcyMjI5MiAxNDguMjIxMDcsMzguMzY3MjIwMyBjIDE0Ny45NjI4MTQsMzguMTA4OTY0NiAxNDguMDQ1ODk4LDM4LjA3NjcyOTkgMTQ3Ljc4NjU2MywzNy45NTQ4Mzc1IGMgMTQ3LjQ3NTUyNywzNy44MDg2NDQxIDE0Ni45Nzc5NDUsMzcuNzEyNTUxNyAxNDYuNjA2NzI5LDM3LjQ2NzcwNzQgYyAxNDYuMjE4Nzk2LDM3LjIxMTgzNzIgMTQ1Ljg0NzE0LDM2LjgxMDY0ODIgMTQ1LjQ2NTE2NywzNi40MDI5NjMyIGMgMTQ1LjE3Njg3MSwzNi4wOTUyNjEzIDE0NC44ODI2OTksMzUuNzgzODU4OSAxNDQuNTcxMjE2LDM1LjUyODQ0MjkgYyAxNDQuMzMwMzg3LDM1LjMzMDk2MjcgMTQzLjczMTE3LDM1LjE2OTA5MjkgMTQzLjMwNDMzNSwzNS4wMzYwMDI3IGMgMTQyLjg3NzUsMzQuOTAyOTEyNSAxNDMuMDgwOSwzNC4wNTc3Nzg0IDE0Mi44NzM5MzYsMzMuODE1NDI3MSBjIDE0Mi42MDA2NjcsMzMuNDk1NDMzOCAxNDIuMTQ4MzU5LDMzLjM3OTAyMzMgMTQxLjQ0NjkxOSwzMy40MTM4ODE1IGMgMTQxLjA0Nzc0LDMzLjQzMzcxODggMTQwLjU2Nzg3OCwzMy41MDI1NDU4IDEzOS45OTQ0MTUsMzMuNjEwNzIwNyBjIDEzOS42MjQ4NywzMy42ODA0Mjk2IDEzOS42MTc4NDMsMzQuMDk0NTcxNCAxMzguOTUxNTk5LDM0LjU4ODQ2OTEgbCAxMzguODExNTE3LDM1LjEwNjEwMTEgeiBcIjtcbi8vIHZhciBwYXRoID0gcG9pbnQxICsgcG9pbnQyICsgcG9pbnQzICsgcG9pbnQ0ICsgcG9pbnQ1ICsgcG9pbnQ2ICsgcG9pbnQ3ICsgcG9pbnQ4ICsgcG9pbnQ5O1xudmFyIHBhdGggPSBsZXR0cmVzQU5UICsgbGV0dHJlUztcblxubW9kdWxlLmV4cG9ydHMgPSBwYXRoOyIsIid1c2Ugc3RyaWN0JztcblxudmFyIHNxcnQgPSBNYXRoLnNxcnQ7XG52YXIgcG93ID0gTWF0aC5wb3c7XG5cbmZ1bmN0aW9uIHNpZ24oeCkge1xuXHRyZXR1cm4geCA/IHggPCAwID8gLTEgOiAxIDogMDtcbn1cblxuZnVuY3Rpb24gcmFuZ2Uoc3RhcnQsIGNvdW50KSB7XG4gICAgcmV0dXJuIEFycmF5LmFwcGx5KDAsIEFycmF5KGNvdW50KSkubWFwKGZ1bmN0aW9uIChlbGVtZW50LCBpbmRleCkge1xuICAgIFx0cmV0dXJuIGluZGV4ICsgc3RhcnRcbiAgICB9KTtcbn1cblxuZnVuY3Rpb24gZGlzdGFuY2UoYSwgYil7XG5cdHJldHVybiBzcXJ0KHBvdyhhLnggLSBiLngsIDIpICsgcG93KGEueSAtIGIueSwgMikpO1xufVxuXG5mdW5jdGlvbiBub3JtKHYpe1xuXHRyZXR1cm4gc3FydChwb3codi54LCAyKSArIHBvdyh2LnksIDIpKTtcbn1cblxubW9kdWxlLmV4cG9ydHMgPSB7XG5cdHNpZ246IHNpZ24sXG5cdHJhbmdlOiByYW5nZSxcblx0ZGlzdGFuY2U6IGRpc3RhbmNlLFxuXHRub3JtOiBub3JtXG59IiwiJ3VzZSBzdHJpY3QnXG5cbmZ1bmN0aW9uIFZlY3Rvcih4LCB5KSB7XG4gICAgdGhpcy54ID0geDsgICAgICAgICAgICAgICAgXG4gICAgdGhpcy55ID0geTtcbn1cblxuVmVjdG9yLnByb3RvdHlwZS5ub3JtID0gZnVuY3Rpb24oKXtcblx0cmV0dXJuIE1hdGguc3FydCh0aGlzLnggKiB0aGlzLnggKyB0aGlzLnkgKiB0aGlzLnkpO1xufVxuXG5WZWN0b3IucHJvdG90eXBlLm5vcm1hbGl6ZSA9IGZ1bmN0aW9uKCl7XG5cdHZhciBub3JtID0gdGhpcy5ub3JtKCk7XG5cdHRoaXMueCA9IHRoaXMueCAvIG5vcm07XG5cdHRoaXMueSA9IHRoaXMueSAvIG5vcm07XG59XG5cblxuXG5tb2R1bGUuZXhwb3J0cyA9IFZlY3RvcjsiLCIndXNlIHN0cmljdCc7XG5cbnZhciBfYW50Q29sb255ID0gcmVxdWlyZSgnLi9pbmRleC5qcycpO1xuXG52YXIgY29udGFpbmVyID0gZG9jdW1lbnQucXVlcnlTZWxlY3RvcignLmNvbG9ueScpO1xuXG52YXIgb3B0aW9ucyA9IHtcblx0dmVsb2NpdHk6IDAuMDAxLFxuXHRuYkFudHM6IDQwMDAsXG5cdHdlaWdodDogMTAsXG5cdHJlcFNpemU6IDAuMDUsXG5cdHJlcFNwZWVkOiAwLjAwMixcblx0bmJTdGFydDogMzAwLFxuXHRuYlJhbmQ6IDMwMFxuXHQvLyBvYmogcGFyIGRlZmF1dFxufTtcblxudmFyIGFudENvbG9ueSA9IF9hbnRDb2xvbnkoY29udGFpbmVyLCBvcHRpb25zKTtcblxud2luZG93LmFkZEV2ZW50TGlzdGVuZXIoJ2NsaWNrJywgZnVuY3Rpb24gKCl7XG5cdC8vIG9wdGlvbnMudmVsb2NpdHkgPSAwLjAwMztcblx0b3B0aW9ucy5uYkFudHMgPSAyMDAwMDtcblx0Ly8gb3B0aW9ucy53ZWlnaHQgPSAxMDAwMDAwMDtcblx0Ly8gb3B0aW9ucy5yZXBTcGVlZCA9IDAuMDE7XG5cdC8vIG9wdGlvbnMucmVwU2l6ZSA9IDAuMTtcblxuXHQvLyBhbnRDb2xvbnkuY2hhbmdlT3B0aW9ucyhvcHRpb25zKTtcblx0YW50Q29sb255LmNoYW5nZU9wdGlvbnMob3B0aW9ucyk7XG59KTtcblxuIl19
