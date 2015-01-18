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
            if (random() > 0.5){
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

        console.log(points);

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

        console.log(points);

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
	var FPSMonitor = document.querySelector('#FPS');
	var dTMonitor = document.querySelector('#dT');
	var refreshTime = 0;
	var maxDeltaTime = 30;
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
			FPSMonitor.style.color = "green";
			population = antsGroup.create(population);
		}	
		else if (antNumber > objPopulation){
			population = antsGroup.remove(population, antNumber - objPopulation);
			FPSMonitor.style.color = "red";
		}
		else
			FPSMonitor.style.color = "white";
	}

	function displayFPS(dT){
		FPSCount = (1000/dT).toFixed(2);
		var t = dT.toFixed(2);
		FPSMonitor.innerText = 'FPS : ' + FPSCount;  
		dTMonitor.innerText = 'nbAnts : ' + population.length;
		// dTMonitor.innerText = 'dT : ' + t + 'ms';
	}

	function tick() {
		var now = performance.now();
		deltaTime = now - lastUpdate;
		lastUpdate = now;
		refreshTime += deltaTime/1000; // in seconds

		console.log('nbAnts', population.length);

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
			objPopulation += 100;
			// FPSOverLimitCount = 0;
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
var lettresANT = "m 25.6185414,31.9447266 l 25.6185414,31.1003906 l 25.388268,30.4095703 l 24.9277211,30.0257813 l 24.5439321,29.5652344 l 24.0833852,29.1814453 l 23.6995961,28.7208984 l 23.0087758,28.490625 l 20.475768,28.490625 l 19.8617055,28.7208984 l 19.4779164,29.1814453 l 18.7870961,29.3349609 l 14.5654164,29.3349609 l 13.9513539,29.5652344 l 13.5675649,30.0257813 l 12.8767446,30.1792969 l 12.0324086,30.1792969 l 11.4183461,30.4095703 l 11.0345571,30.8701172 l 10.3437367,31.0236328 l 9.72967424,31.2539063 l 9.34588518,31.7144531 l 8.8853383,32.0982422 l 8.50154924,32.5587891 l 8.04100237,32.9425781 l 7.88748674,33.6333984 l 7.6572133,34.2474609 l 7.19666643,34.63125 l 7.0431508,35.3220703 l 7.0431508,36.1664063 l 7.19666643,36.7804688 l 7.6572133,37.1642578 l 8.04100237,37.6248047 l 8.50154924,38.0085938 l 8.8853383,38.4691406 l 9.34588518,38.8529297 l 9.72967424,39.3134766 l 10.3437367,39.4669922 l 15.4097524,39.4669922 l 16.1005727,39.3134766 l 16.4843617,38.8529297 l 17.0984242,38.6226563 l 17.7892446,38.4691406 l 18.1730336,38.0085938 l 18.6335805,37.6248047 l 19.0173696,37.1642578 l 19.6314321,36.9339844 l 20.3222524,36.7804688 l 20.7060414,36.3199219 l 21.3201039,36.0896484 l 22.0109242,35.9361328 l 22.3947133,35.4755859 l 23.0087758,35.2453125 l 23.6995961,35.0917969 l 23.9298696,34.4777344 l 24.0833852,33.7869141 l 24.5439321,33.403125 l 24.9277211,32.9425781 l 25.388268,32.5587891 l 25.6185414,31.9447266  m 37.5927602,43.9189453 l 37.2089711,44.3794922 l 36.7484242,44.7632813 l 36.3646352,45.2238281 l 35.6738149,45.3773438 l 35.0597524,45.6076172 l 34.6759633,46.0681641 l 33.985143,46.2216797 l 31.4521352,46.2216797 l 30.8380727,46.0681641 l 30.4542836,45.6076172 l 29.7634633,45.3773438 l 29.1494008,45.2238281 l 28.7656117,44.7632813 l 28.3050649,44.3794922 l 27.9212758,43.9189453 l 27.2304555,43.6886719 l 26.616393,43.5351563 l 26.2326039,43.0746094 l 25.7720571,43.0746094 l 25.388268,43.5351563 l 24.6974477,43.6886719 l 23.8531117,43.6886719 l 23.2390492,43.9189453 l 22.8552602,44.3794922 l 22.1644399,44.5330078 l 21.5503774,44.7632813 l 21.1665883,45.2238281 l 20.475768,45.3773438 l 19.6314321,45.3773438 l 19.0173696,45.6076172 l 18.6335805,46.0681641 l 17.9427602,46.2216797 l 16.2540883,46.2216797 l 15.6400258,46.4519531 l 15.2562367,46.9125 l 14.7956899,46.9125 l 14.4119008,46.4519531 l 13.7210805,46.2216797 l 9.4994008,46.2216797 l 8.8853383,46.0681641 l 8.50154924,45.6076172 l 7.81072893,45.3773438 l 7.19666643,45.2238281 l 6.81287737,44.7632813 l 6.12205705,44.5330078 l 5.50799455,44.3794922 l 5.12420549,43.9189453 l 4.66365862,43.5351563 l 4.27986955,43.0746094 l 3.81932268,42.6908203 l 3.43553362,42.2302734 l 2.97498674,41.8464844 l 2.59119768,41.3859375 l 2.1306508,41.0021484 l 1.97713518,40.3880859 l 1.74686174,39.6972656 l 1.28631487,39.3134766 l 1.13279924,38.6994141 l 1.13279924,37.8550781 l 0.902525804,37.1642578 l 0.441978929,36.7804688 l 0.288463304,36.1664063 l 0.288463304,33.6333984 l 0.441978929,32.9425781 l 0.902525804,32.5587891 l 1.13279924,31.9447266 l 1.13279924,31.1003906 l 1.28631487,30.4095703 l 1.74686174,30.0257813 l 1.97713518,29.4117188 l 2.1306508,28.7208984 l 2.59119768,28.3371094 l 2.97498674,27.8765625 l 3.43553362,27.4927734 l 3.81932268,27.0322266 l 4.27986955,26.6484375 l 4.66365862,26.1878906 l 5.27772112,25.9576172 l 6.12205705,25.9576172 l 6.81287737,25.8041016 l 7.19666643,25.3435547 l 7.81072893,25.1132813 l 8.65506487,25.1132813 l 9.34588518,24.9597656 l 9.72967424,24.4992188 l 10.3437367,24.2689453 l 13.7210805,24.2689453 l 14.4119008,24.1154297 l 14.7956899,23.6548828 l 15.4097524,23.4246094 l 17.9427602,23.4246094 l 18.6335805,23.2710938 l 19.0173696,22.8105469 l 19.6314321,22.5802734 l 23.8531117,22.5802734 l 24.5439321,22.4267578 l 24.9277211,21.9662109 l 25.388268,21.5824219 l 25.7720571,21.121875 l 26.3861196,20.8916016 l 27.0769399,20.7380859 l 27.3072133,20.1240234 l 27.3072133,18.4353516 l 27.0769399,17.7445313 l 26.616393,17.3607422 l 26.616393,16.9001953 l 27.0769399,16.5164063 l 27.0769399,16.0558594 l 26.616393,15.6720703 l 26.4628774,15.0580078 l 26.2326039,14.3671875 l 25.7720571,13.9833984 l 25.388268,13.5228516 l 24.9277211,13.1390625 l 24.5439321,12.6785156 l 23.8531117,12.4482422 l 23.0087758,12.4482422 l 22.3947133,12.2947266 l 22.0109242,11.8341797 l 21.5503774,11.8341797 l 21.1665883,12.2947266 l 20.475768,12.4482422 l 14.5654164,12.4482422 l 13.9513539,12.6785156 l 13.5675649,13.1390625 l 13.107018,13.5228516 l 12.7232289,13.9833984 l 12.2626821,14.3671875 l 12.1091664,15.0580078 l 12.1091664,17.5910156 l 11.878893,18.2050781 l 11.4183461,18.5888672 l 11.0345571,19.0494141 l 10.5740102,19.4332031 l 10.1902211,19.89375 l 9.4994008,20.0472656 l 7.81072893,20.0472656 l 7.19666643,19.89375 l 6.81287737,19.4332031 l 6.12205705,19.2029297 l 5.50799455,19.0494141 l 5.12420549,18.5888672 l 4.66365862,18.2050781 l 4.51014299,17.5910156 l 4.51014299,16.7466797 l 4.27986955,16.0558594 l 3.81932268,15.6720703 l 3.66580705,15.0580078 l 3.66580705,14.2136719 l 3.81932268,13.5228516 l 4.27986955,13.1390625 l 4.51014299,12.525 l 4.66365862,11.8341797 l 5.12420549,11.4503906 l 5.50799455,10.9898438 l 5.96854143,10.6060547 l 6.35233049,10.1455078 l 6.96639299,9.91523438 l 7.6572133,9.76171875 l 8.04100237,9.30117188 l 8.65506487,9.07089844 l 9.34588518,8.91738281 l 9.72967424,8.45683594 l 10.3437367,8.2265625 l 13.7210805,8.2265625 l 14.4119008,8.07304688 l 14.7956899,7.6125 l 15.2562367,7.6125 l 15.6400258,8.07304688 l 16.2540883,8.2265625 l 19.6314321,8.2265625 l 20.3222524,8.45683594 l 20.7060414,8.91738281 l 21.3201039,9.07089844 l 23.8531117,9.07089844 l 24.5439321,9.30117188 l 24.9277211,9.76171875 l 25.5417836,9.91523438 l 26.2326039,10.1455078 l 26.616393,10.6060547 l 27.0769399,10.9898438 l 27.4607289,11.4503906 l 27.9212758,11.8341797 l 28.3050649,12.2947266 l 28.7656117,12.6785156 l 29.1494008,13.1390625 l 29.6099477,13.5228516 l 29.8402211,14.2136719 l 29.9937367,14.8277344 l 30.4542836,15.2115234 l 30.6845571,15.9023438 l 30.6845571,20.1240234 l 30.8380727,20.7380859 l 31.2986196,21.121875 l 31.528893,21.8126953 l 31.528893,35.3220703 l 31.6824086,35.9361328 l 32.1429555,36.3199219 l 32.5267446,36.7804688 l 32.9872914,37.1642578 l 33.3710805,37.6248047 l 33.985143,37.7783203 l 34.6759633,38.0085938 l 35.0597524,38.4691406 l 35.6738149,38.6226563 l 36.3646352,38.8529297 l 36.7484242,39.3134766 l 37.3624867,39.4669922 l 38.0533071,39.6972656 l 38.2835805,40.3880859 l 38.2835805,42.9210938 l 38.0533071,43.5351563 l 37.5927602,43.9189453 z m 87.4181857,40.5416016 l 87.2646701,41.2324219 l 87.0343966,41.8464844 l 86.5738497,42.2302734 l 86.1900607,42.6908203 l 85.7295138,43.0746094 l 85.3457247,43.5351563 l 84.6549044,43.6886719 l 84.0408419,43.9189453 l 83.6570529,44.3794922 l 82.9662326,44.5330078 l 72.8342013,44.5330078 l 72.2201388,44.3794922 l 71.8363497,43.9189453 l 71.3758029,43.5351563 l 71.2222872,42.9210938 l 70.9920138,42.2302734 l 70.5314669,41.8464844 l 70.3779513,41.2324219 l 70.5314669,40.5416016 l 70.9920138,40.1578125 l 71.3758029,39.6972656 l 71.8363497,39.3134766 l 72.2201388,38.8529297 l 72.6806857,38.4691406 l 73.0644747,38.0085938 l 73.6785372,37.7783203 l 74.3693576,37.6248047 l 74.7531466,37.1642578 l 75.2136935,36.7804688 l 75.5974826,36.3199219 l 76.0580294,35.9361328 l 76.2883029,35.3220703 l 76.2883029,34.4777344 l 76.4418185,33.7869141 l 76.9023654,33.403125 l 77.1326388,32.7890625 l 77.1326388,17.5910156 l 77.2861544,16.9001953 l 77.7467013,16.5164063 l 77.9769747,15.9023438 l 77.7467013,15.2115234 l 77.2861544,14.8277344 l 77.1326388,14.2136719 l 76.9023654,13.5228516 l 76.4418185,13.1390625 l 76.0580294,12.6785156 l 75.3672091,12.4482422 l 74.7531466,12.2947266 l 74.3693576,11.8341797 l 73.6785372,11.6039063 l 73.0644747,11.4503906 l 72.6806857,10.9898438 l 71.9898654,10.7595703 l 66.0795138,10.7595703 l 65.4654513,10.9898438 l 65.0816622,11.4503906 l 64.3908419,11.6039063 l 63.7767794,11.8341797 l 63.6232638,12.525 l 63.3929904,13.1390625 l 62.9324435,13.5228516 l 62.5486544,13.9833984 l 62.0881076,14.3671875 l 61.7043185,14.8277344 l 61.2437716,15.2115234 l 61.090256,15.9023438 l 60.8599826,16.5164063 l 60.3994357,16.9001953 l 60.0156466,17.3607422 l 59.5550997,17.7445313 l 59.4015841,18.4353516 l 59.1713107,19.0494141 l 58.7107638,19.4332031 l 58.5572482,20.1240234 l 58.3269747,20.7380859 l 57.8664279,21.121875 l 57.7129122,21.8126953 l 57.4826388,22.4267578 l 57.0220919,22.8105469 l 56.8685763,23.5013672 l 56.8685763,28.5673828 l 57.0220919,29.1814453 l 57.4826388,29.5652344 l 57.7129122,30.2560547 l 57.7129122,31.1003906 l 57.4826388,31.7144531 l 57.0220919,32.0982422 l 56.8685763,32.7890625 l 57.0220919,33.403125 l 57.4826388,33.7869141 l 57.7129122,34.4777344 l 57.7129122,35.3220703 l 57.8664279,35.9361328 l 58.3269747,36.3199219 l 58.7107638,36.7804688 l 59.1713107,37.1642578 l 59.5550997,37.6248047 l 60.1691622,37.7783203 l 60.8599826,38.0085938 l 61.2437716,38.4691406 l 61.7043185,38.8529297 l 62.0881076,39.3134766 l 62.5486544,39.6972656 l 62.7789279,40.3880859 l 62.9324435,41.0021484 l 63.3929904,41.3859375 l 63.3929904,41.8464844 l 62.9324435,42.2302734 l 62.7789279,42.9210938 l 62.5486544,43.5351563 l 61.8578341,43.6886719 l 61.2437716,43.9189453 l 60.8599826,44.3794922 l 60.1691622,44.5330078 l 59.3248263,44.5330078 l 58.7107638,44.7632813 l 58.3269747,45.2238281 l 57.8664279,45.2238281 l 57.4826388,44.7632813 l 57.0220919,44.7632813 l 56.6383029,45.2238281 l 55.9474826,45.3773438 l 55.1031466,45.3773438 l 54.4890841,45.2238281 l 54.1052951,44.7632813 l 53.6447482,44.7632813 l 53.2609591,45.2238281 l 52.5701388,45.3773438 l 49.1927951,45.3773438 l 48.5787326,45.2238281 l 48.1949435,44.7632813 l 47.5041232,44.5330078 l 46.8900607,44.3794922 l 46.5062716,43.9189453 l 46.0457247,43.5351563 l 45.6619357,43.0746094 l 45.2013888,42.6908203 l 45.0478732,42.0767578 l 45.0478732,41.2324219 l 45.2013888,40.5416016 l 45.6619357,40.1578125 l 46.0457247,39.6972656 l 46.5062716,39.3134766 l 46.8900607,38.8529297 l 47.5041232,38.6226563 l 48.1949435,38.4691406 l 48.5787326,38.0085938 l 49.0392794,37.6248047 l 49.4230685,37.1642578 l 49.8836154,36.7804688 l 50.2674044,36.3199219 l 50.7279513,35.9361328 l 50.9582247,35.3220703 l 50.9582247,26.034375 l 51.1117404,25.3435547 l 51.5722872,24.9597656 l 51.8025607,24.3457031 l 51.8025607,23.5013672 l 51.5722872,22.8105469 l 51.1117404,22.4267578 l 50.9582247,21.8126953 l 50.9582247,16.7466797 l 50.7279513,16.0558594 l 50.2674044,15.6720703 l 50.1138888,15.0580078 l 49.8836154,14.3671875 l 49.4230685,13.9833984 l 49.0392794,13.5228516 l 48.5787326,13.1390625 l 48.1949435,12.6785156 l 47.7343966,12.2947266 l 47.3506076,11.8341797 l 46.8900607,11.4503906 l 46.5062716,10.9898438 l 45.8154513,10.7595703 l 45.2013888,10.6060547 l 45.2013888,10.1455078 l 45.6619357,9.76171875 l 45.8922091,9.14765625 l 45.8922091,8.30332031 l 46.0457247,7.6125 l 46.5062716,7.22871094 l 46.8900607,6.76816406 l 47.3506076,6.384375 l 47.7343966,5.92382813 l 48.3484591,5.69355469 l 49.0392794,5.54003906 l 49.4230685,5.07949219 l 50.037131,4.84921875 l 50.7279513,5.07949219 l 51.1117404,5.54003906 l 51.7258029,5.69355469 l 52.5701388,5.69355469 l 53.2609591,5.92382813 l 53.6447482,6.384375 l 54.2588107,6.53789063 l 54.949631,6.76816406 l 55.3334201,7.22871094 l 55.7939669,7.6125 l 56.177756,8.07304688 l 56.6383029,8.45683594 l 57.0220919,8.91738281 l 57.4826388,9.30117188 l 57.8664279,9.76171875 l 58.3269747,10.1455078 l 58.7107638,10.6060547 l 59.1713107,10.6060547 l 59.5550997,10.1455078 l 60.1691622,9.91523438 l 61.0134982,9.91523438 l 61.7043185,9.76171875 l 62.0881076,9.30117188 l 62.5486544,8.91738281 l 62.9324435,8.45683594 l 63.3929904,8.07304688 l 63.7767794,7.6125 l 64.3908419,7.38222656 l 65.0816622,7.22871094 l 65.4654513,6.76816406 l 66.0795138,6.53789063 l 68.6125216,6.53789063 l 69.3033419,6.384375 l 69.687131,5.92382813 l 70.3011935,5.69355469 l 73.6785372,5.69355469 l 74.3693576,5.92382813 l 74.7531466,6.384375 l 75.3672091,6.53789063 l 76.0580294,6.76816406 l 76.4418185,7.22871094 l 77.055881,7.38222656 l 77.7467013,7.6125 l 78.1304904,8.07304688 l 78.5910372,8.45683594 l 78.9748263,8.91738281 l 79.4353732,9.30117188 l 79.8191622,9.76171875 l 80.2797091,10.1455078 l 80.5099826,10.8363281 l 80.6634982,11.4503906 l 81.1240451,11.8341797 l 81.3543185,12.525 l 81.5078341,13.1390625 l 81.968381,13.5228516 l 82.1986544,14.2136719 l 82.1986544,16.7466797 l 81.968381,17.3607422 l 81.5078341,17.7445313 l 81.5078341,18.2050781 l 81.968381,18.5888672 l 81.968381,19.0494141 l 81.5078341,19.4332031 l 81.3543185,20.1240234 l 81.3543185,21.8126953 l 81.5078341,22.4267578 l 81.968381,22.8105469 l 82.1986544,23.5013672 l 82.1986544,28.5673828 l 82.3521701,29.1814453 l 82.8127169,29.5652344 l 83.0429904,30.2560547 l 83.0429904,31.1003906 l 82.8127169,31.7144531 l 82.3521701,32.0982422 l 82.1986544,32.7890625 l 82.3521701,33.403125 l 82.8127169,33.7869141 l 83.0429904,34.4777344 l 83.196506,35.0917969 l 83.6570529,35.4755859 l 84.0408419,35.9361328 l 84.6549044,36.0896484 l 85.3457247,36.3199219 l 85.7295138,36.7804688 l 86.3435763,36.9339844 l 87.0343966,37.1642578 l 87.2646701,37.8550781 l 87.4181857,38.4691406 l 87.8787326,38.8529297 l 88.109006,39.54375 l 87.8787326,40.1578125 l 87.4181857,40.5416016 z m 130.488924,32.0982422 l 130.335408,32.7890625 l 130.335408,34.4777344 l 130.105135,35.0917969 l 129.644588,35.4755859 l 129.491072,36.1664063 l 129.491072,37.8550781 l 129.260799,38.4691406 l 128.800252,38.8529297 l 128.646736,39.54375 l 128.416463,40.1578125 l 127.955916,40.5416016 l 127.572127,41.0021484 l 127.11158,41.3859375 l 126.727791,41.8464844 l 126.267244,42.2302734 l 125.883455,42.6908203 l 125.422908,43.0746094 l 125.039119,43.5351563 l 124.578572,43.9189453 l 124.194783,44.3794922 l 123.503963,44.5330078 l 122.8899,44.7632813 l 122.506111,45.2238281 l 121.815291,45.3773438 l 120.970955,45.3773438 l 120.356892,45.6076172 l 119.973103,46.0681641 l 119.282283,46.2216797 l 115.060603,46.2216797 l 114.446541,46.0681641 l 114.062752,45.6076172 l 113.371932,45.3773438 l 112.527596,45.3773438 l 111.913533,45.2238281 l 111.529744,44.7632813 l 110.838924,44.5330078 l 110.224861,44.3794922 l 109.841072,43.9189453 l 109.380525,43.5351563 l 108.996736,43.0746094 l 108.536189,42.6908203 l 108.1524,42.2302734 l 107.691853,41.8464844 l 107.308064,41.3859375 l 106.847517,41.0021484 l 106.463728,40.5416016 l 106.003182,40.1578125 l 105.619392,39.6972656 l 105.158846,39.3134766 l 105.00533,38.6994141 l 104.775057,38.0085938 l 104.31451,37.6248047 l 104.160994,37.0107422 l 104.160994,36.1664063 l 103.930721,35.4755859 l 103.470174,35.0917969 l 103.316658,34.4777344 l 103.316658,32.7890625 l 103.470174,32.0982422 l 103.930721,31.7144531 l 104.160994,31.1003906 l 104.160994,27.7230469 l 104.31451,27.0322266 l 104.775057,26.6484375 l 105.00533,26.034375 l 104.775057,25.3435547 l 104.31451,24.9597656 l 104.160994,24.3457031 l 104.31451,23.6548828 l 104.775057,23.2710938 l 105.00533,22.6570313 l 105.00533,15.9023438 l 104.775057,15.2115234 l 104.31451,14.8277344 l 104.160994,14.2136719 l 103.930721,13.5228516 l 103.470174,13.1390625 l 103.086385,12.6785156 l 102.395564,12.4482422 l 101.781502,12.2947266 l 101.397713,11.8341797 l 100.706892,11.6039063 l 98.1738846,11.6039063 l 97.5598221,11.4503906 l 97.1760331,10.9898438 l 96.7154862,10.6060547 l 96.3316971,10.1455078 l 95.8711503,9.76171875 l 95.4873612,9.30117188 l 95.0268143,8.91738281 l 94.8732987,8.30332031 l 95.0268143,7.6125 l 95.4873612,7.22871094 l 95.8711503,6.76816406 l 96.3316971,6.384375 l 96.7154862,5.92382813 l 97.1760331,5.54003906 l 97.5598221,5.07949219 l 98.1738846,4.84921875 l 102.395564,4.84921875 l 103.086385,5.07949219 l 103.470174,5.54003906 l 103.930721,5.54003906 l 104.31451,5.07949219 l 104.775057,4.69570313 l 105.00533,4.08164063 l 105.158846,3.39082031 l 105.619392,3.00703125 l 105.849666,2.39296875 l 106.003182,1.70214844 l 106.617244,1.471875 l 107.308064,1.31835938 l 107.538338,0.704296875 l 107.691853,0.0134765625 l 108.1524,0.0134765625 l 108.382674,0.704296875 l 108.536189,1.31835938 l 108.996736,1.70214844 l 109.380525,2.16269531 l 109.841072,2.54648438 l 110.071346,3.23730469 l 110.071346,4.08164063 l 110.224861,4.69570313 l 110.685408,5.07949219 l 111.069197,5.54003906 l 111.68326,5.69355469 l 112.37408,5.92382813 l 112.757869,6.384375 l 113.218416,6.384375 l 113.602205,5.92382813 l 114.216267,5.69355469 l 115.904939,5.69355469 l 116.59576,5.92382813 l 116.979549,6.384375 l 117.593611,6.53789063 l 122.659627,6.53789063 l 123.350447,6.76816406 l 123.734236,7.22871094 l 124.194783,7.6125 l 124.578572,8.07304688 l 125.039119,8.45683594 l 125.269392,9.14765625 l 125.039119,9.76171875 l 124.578572,10.1455078 l 124.194783,10.6060547 l 123.734236,10.9898438 l 123.350447,11.4503906 l 122.659627,11.6039063 l 113.371932,11.6039063 l 112.757869,11.8341797 l 112.37408,12.2947266 l 111.913533,12.6785156 l 111.529744,13.1390625 l 111.069197,13.5228516 l 110.685408,13.9833984 l 110.224861,14.3671875 l 110.071346,15.0580078 l 110.071346,16.7466797 l 110.224861,17.3607422 l 110.685408,17.7445313 l 110.685408,18.2050781 l 110.224861,18.5888672 l 110.071346,19.2796875 l 110.071346,28.5673828 l 110.224861,29.1814453 l 110.685408,29.5652344 l 110.915682,30.2560547 l 110.915682,31.1003906 l 111.069197,31.7144531 l 111.529744,32.0982422 l 111.760017,32.7890625 l 111.760017,33.6333984 l 111.913533,34.2474609 l 112.37408,34.63125 l 112.604353,35.3220703 l 112.757869,35.9361328 l 113.371932,36.0896484 l 114.062752,36.3199219 l 114.446541,36.7804688 l 115.060603,36.9339844 l 118.437947,36.9339844 l 119.128767,36.7804688 l 119.512557,36.3199219 l 119.973103,35.9361328 l 120.356892,35.4755859 l 120.817439,35.0917969 l 121.047713,34.4777344 l 121.201228,33.7869141 l 121.661775,33.403125 l 121.892049,32.7890625 l 121.892049,31.9447266 l 122.045564,31.2539063 l 122.506111,30.8701172 l 122.736385,30.2560547 l 122.8899,29.5652344 l 123.350447,29.1814453 l 123.734236,28.7208984 l 124.348299,28.490625 l 125.039119,28.3371094 l 125.422908,27.8765625 l 126.036971,27.6462891 l 126.881307,27.6462891 l 127.572127,27.8765625 l 127.955916,28.3371094 l 128.569978,28.490625 l 129.260799,28.7208984 l 129.644588,29.1814453 l 130.105135,29.5652344 l 130.335408,30.2560547 l 130.488924,30.8701172 l 130.949471,31.2539063 l 130.949471,31.7144531 l 130.488924,32.0982422 z";
var lettreS = "m 138.811517,35.1061011 c 138.63781,35.7365111 138.606564,35.8010614 138.198455,36.2091704 c 138.056753,36.3508722 137.999155,37.0215928 137.948798,37.2230199 c 137.869537,37.5400622 137.926411,37.8292331 137.948798,38.0775475 c 137.983078,38.4577806 137.927766,38.9415367 138.126037,39.2217479 c 138.240713,39.3838164 138.438596,39.433835 138.560544,39.6201964 c 138.747066,39.9052408 138.79591,40.4584527 138.951599,40.833933 c 139.061131,41.0980922 139.263918,41.2045958 139.429556,41.3702343 c 139.860819,41.8014965 139.387712,43.8969432 139.699324,44.8989806 c 139.798518,45.2179557 140.058533,45.3368515 140.249382,45.5162438 c 140.468372,45.7220877 140.667698,45.9527405 140.883714,46.2064196 c 141.287987,46.6811791 141.363563,46.8042662 141.758733,46.8889452 c 142.153903,46.9736241 142.389571,46.8889452 142.5773,46.8889452 c 142.919018,46.8889452 143.062586,46.5447924 143.170055,46.4373229 c 143.5244,46.0829776 144.30261,46.0849327 144.456107,45.9576928 c 144.609603,45.8304528 144.538252,45.788898 144.680169,45.6469811 c 145.235946,45.0912041 148.31368,45.2113478 149.221661,45.3929442 c 149.793527,45.5073172 149.592293,45.9406794 150.311344,46.0844896 c 150.836944,46.1896096 151.765006,46.0703795 152.400101,46.2064196 c 153.035195,46.3424596 152.843259,46.6989649 153.303797,46.8553563 c 154.241306,47.173719 154.949065,47.0720881 155.420329,46.8042662 c 155.701153,46.6446723 155.856582,46.3018852 156.238444,46.2064196 c 157.795519,45.8171509 159.745397,46.501263 161.376107,45.9576928 c 161.556165,45.8976734 161.595501,45.6542159 161.855957,45.5162438 c 162.214357,45.3263872 162.76045,45.2045722 162.999127,45.1250133 c 163.301469,45.0242326 163.254211,44.8176843 163.542798,44.673391 c 163.811921,44.5388293 164.198751,44.5040294 164.537469,44.3628991 c 164.876187,44.2217687 165.245628,43.7958718 165.59596,43.4455407 c 165.694074,43.3474263 166.089427,42.951577 166.414031,42.5846337 c 166.738635,42.2176903 166.649398,41.6958381 166.82381,41.3702343 c 166.942587,41.1484928 167.231728,41.0717202 167.346164,40.833933 c 167.514797,40.4835286 167.470564,39.9584339 167.642542,39.6201964 c 167.752332,39.4042681 168.064756,39.3304646 168.164732,39.1262328 c 168.337576,38.7731432 168.382235,38.2192425 168.517395,37.9548375 c 168.663802,37.6684294 168.946353,37.5366951 169.011358,37.4066857 c 169.322207,36.7849881 169.138191,35.2241166 169.350074,34.5884691 c 169.469774,34.2293693 169.81581,34.2906902 169.971058,33.7685027 c 170.126306,33.2463152 169.678613,32.7676354 169.40929,32.4983129 c 169.209305,32.2983273 169.341035,31.4437653 168.57385,30.8470659 c 168.319811,30.6494805 168.456555,30.0086667 168.164732,29.6053234 c 167.993437,29.3685691 167.724169,29.2856739 167.642542,29.0407946 c 167.455896,28.4808568 167.609275,28.4279622 167.346164,27.982304 c 167.124194,27.6063293 166.62278,27.1459759 166.287673,26.8554088 c 165.944635,26.5579643 165.711054,26.2484947 165.454331,26.0064531 c 165.144628,25.7144613 165.045061,25.484594 164.678766,25.2711034 c 164.34227,25.0749805 163.922398,25.0466617 163.542798,24.8773955 c 163.194478,24.7220773 163.083972,24.3552062 162.507536,24.2281877 c 161.9311,24.1011693 161.486505,24.2278936 161.023442,24.0729418 c 160.683926,23.9593313 160.613722,23.5600169 160.247214,23.4378478 c 159.14903,23.0717863 156.987132,23.713421 155.927908,23.1838091 c 155.706656,23.0731829 155.491061,22.6254308 154.912418,22.5487133 c 154.333776,22.4719957 153.815269,22.5408899 153.477743,22.3953734 c 153.162857,22.2596177 153.059763,21.8669612 152.90768,21.8289403 c 150.954444,21.3406314 149.563466,21.9550928 148.358942,21.4371413 c 147.975987,21.272469 147.898826,21.0950527 147.20865,20.4176183 c 146.518475,19.7401839 145.985668,19.3309002 145.738826,18.9639568 c 145.491984,18.5970134 145.583581,18.1453902 145.368182,17.7925597 c 145.152783,17.4397292 145.165926,17.7984953 144.896094,16.9889995 c 144.810593,16.7324952 145.287599,16.6014434 145.368182,16.2988237 c 145.639372,15.2804067 145.643957,15.1850966 146.165982,14.6594246 c 146.360127,14.4639228 146.490629,13.4976765 146.72609,13.2622164 c 147.187271,12.8010353 147.770467,12.442705 148.024506,11.7803283 c 148.278545,11.1179515 147.982167,11.0791053 148.358942,10.9053089 c 148.735717,10.7315125 148.771309,11.160942 149.221661,11.3992716 c 149.504567,11.5489873 149.822314,11.5203687 150.197941,11.5545166 c 151.074442,11.6341985 153.434038,11.4608384 153.61334,11.6401409 c 153.766508,11.7933085 153.971577,12.0864654 154.206096,12.2037249 c 155.25581,12.7285818 156.833988,12.0287261 157.86142,12.5424421 c 157.984728,12.6040962 158.233373,12.9583493 158.341468,13.0123969 c 159.024284,13.3538048 159.278204,13.2615788 159.639685,13.4588758 c 159.942308,13.6240473 160.113559,13.8965364 160.486478,14.2710927 c 160.686801,14.4722946 160.886652,14.685981 161.1847,14.9840285 c 161.61246,15.4117885 162.871104,16.9047033 163.435633,17.2575338 c 164.000162,17.6103643 164.021282,17.325205 164.339387,17.4115026 c 165.215769,17.6492529 166.033115,17.1057564 166.509697,16.2988237 c 166.869078,15.6903318 166.693209,14.5237991 166.682182,13.2622164 c 166.674733,12.4101235 166.775548,11.3232784 166.414031,10.9617617 c 166.128028,10.6757592 166.159076,10.7905709 165.962407,10.3972329 c 165.856212,10.1848444 165.924877,9.59242206 165.454331,9.16938297 c 164.983785,8.74634388 164.974481,8.36180797 164.452292,8.26613685 c 162.570268,7.92132827 160.804751,8.15650513 158.849345,8.15323116 c 158.30471,8.15231927 157.5235,8.19792475 157.02874,8.07423488 c 156.478673,7.93671794 156.410757,7.53795352 155.927908,7.3770041 c 155.448624,7.21724255 154.659909,7.3770041 153.7078,7.27887313 c 153.241169,7.23077885 153.177424,6.91365194 152.90768,6.7284576 c 152.637688,6.54309314 152.23055,6.41794086 151.566924,6.47441966 c 151.12078,6.51238933 150.83207,6.57321218 150.371953,7.12362771 c 149.802771,7.36355241 149.279692,7.3710543 149.012431,7.50468461 c 148.716054,8.49260991 146.88326,7.997023 146.518475,8.36180796 c 146.394821,8.48546169 146.314956,8.65382226 146.165982,8.80279659 c 145.948678,9.02010087 145.440191,8.89079158 145.130308,9.04573276 c 144.88043,9.17067194 144.66026,9.44863746 144.456107,9.65279072 c 144.346843,9.76205412 143.798172,9.67084492 143.304335,10.0297316 c 142.810498,10.3886184 142.157865,11.0432512 142.08985,11.1355582 c 141.957628,11.3150011 141.643424,11.5108743 141.5439,11.7099224 c 141.102752,12.5922176 141.357822,12.4703809 141.236493,12.6980144 c 141.176405,12.8107505 141.061358,13.0475396 140.650022,13.4588759 c 140.642454,13.4664439 140.475326,13.8671445 140.479441,14.2710924 c 140.483232,14.6433395 140.496244,15.0055832 140.409424,15.1792235 c 140.115401,15.7672686 140.079292,15.5288997 139.848007,15.9935344 c 139.757941,16.1744701 139.58036,16.6889278 139.699324,16.9889992 c 139.885854,17.4594979 140.374476,17.9028886 140.409424,18.0776274 c 140.47668,18.413908 140.479334,18.662833 140.479441,19.1170398 c 140.479688,20.1758834 141.4426,20.2917516 141.906073,21.0593792 c 142.369547,21.8270067 142.205884,22.3311986 142.369547,22.3953734 c 142.760603,22.5487133 143.328866,22.4966556 143.600649,22.768438 c 143.883821,23.05161 143.832385,23.9995403 144.093089,24.0729418 c 144.78673,24.2682369 144.91845,24.3326604 145.266426,24.4485295 c 145.614403,24.5643987 145.688025,24.9264879 146.165982,25.3320262 c 146.643938,25.7375646 146.888584,26.0064531 147.20865,26.3748412 c 147.681822,26.9194498 148.443557,26.5631285 148.945247,27.3017888 c 149.353529,27.9029185 151.889664,27.3111578 152.71097,27.7218108 c 153.119776,27.9262138 152.94786,27.9045245 153.190169,28.1468334 c 153.790721,28.7473856 155.332337,28.185207 156.085636,28.5618566 c 156.437966,28.7380211 156.369805,28.9284363 156.940166,29.2136168 c 157.164927,29.3259974 158.200234,29.2136168 158.634741,29.4308696 c 159.069248,29.6481224 159.03255,29.8598897 159.195756,29.9142919 c 159.579797,30.0423054 160.386616,30.1669654 160.609554,30.3899032 c 160.99896,30.7793095 161.359188,31.2508628 161.647319,31.5165013 c 162.189114,32.0160027 162.438162,32.3072349 162.438162,33.6107207 c 162.438162,34.0186619 162.452837,34.7374085 162.226661,34.9635845 c 162.099415,35.0908296 161.34304,35.1740848 161.184699,35.2532555 c 160.803153,35.4440287 160.436636,35.8990815 160.126547,36.2091704 c 159.841794,36.4939233 159.179177,37.2442251 158.849344,37.5271727 c 158.562323,37.7733947 157.94583,37.6887831 157.649859,37.8458129 c 157.288559,38.0375035 157.304334,38.2294232 157.028739,38.3672203 c 156.288135,38.7375227 153.945628,38.3672203 153.477742,38.7003405 c 152.907183,39.1065606 153.023548,39.4969371 151.233648,39.3231334 c 150.418741,39.2440039 150.596106,38.8596605 150.197941,38.7003405 c 149.799777,38.5410204 148.403468,38.5496183 148.22107,38.3672203 c 147.962814,38.1089646 148.045898,38.0767299 147.786563,37.9548375 c 147.475527,37.8086441 146.977945,37.7125517 146.606729,37.4677074 c 145.926002,37.0187181 145.295393,36.1222686 144.571216,35.5284429 c 144.330387,35.3309627 143.73117,35.1690929 143.304335,35.0360027 c 142.8775,34.9029125 143.0809,34.0577784 142.873936,33.8154271 c 142.445154,33.3133302 141.575571,33.3124602 139.994415,33.6107207 c 139.62487,33.6804296 139.617843,34.0945714 138.951599,34.5884691 l 138.811517,35.1061011 z";
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
//# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJzb3VyY2VzIjpbIi91c3IvbG9jYWwvbGliL25vZGVfbW9kdWxlcy93YXRjaGlmeS9ub2RlX21vZHVsZXMvYnJvd3NlcmlmeS9ub2RlX21vZHVsZXMvYnJvd3Nlci1wYWNrL19wcmVsdWRlLmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9hbnRzL0FudENvbG9ueS9pbmRleC5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vYW50cy9BbnRDb2xvbnkvbm9kZV9tb2R1bGVzL2RlbGF1bmF5LXRyaWFuZ3VsYXRlL25vZGVfbW9kdWxlcy9pbmNyZW1lbnRhbC1jb252ZXgtaHVsbC9pY2guanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL2FudHMvQW50Q29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvaW5jcmVtZW50YWwtY29udmV4LWh1bGwvbm9kZV9tb2R1bGVzL3JvYnVzdC1vcmllbnRhdGlvbi9ub2RlX21vZHVsZXMvcm9idXN0LXNjYWxlL25vZGVfbW9kdWxlcy90d28tc3VtL3R3by1zdW0uanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL2FudHMvQW50Q29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvaW5jcmVtZW50YWwtY29udmV4LWh1bGwvbm9kZV9tb2R1bGVzL3JvYnVzdC1vcmllbnRhdGlvbi9ub2RlX21vZHVsZXMvcm9idXN0LXNjYWxlL3JvYnVzdC1zY2FsZS5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vYW50cy9BbnRDb2xvbnkvbm9kZV9tb2R1bGVzL2RlbGF1bmF5LXRyaWFuZ3VsYXRlL25vZGVfbW9kdWxlcy9pbmNyZW1lbnRhbC1jb252ZXgtaHVsbC9ub2RlX21vZHVsZXMvcm9idXN0LW9yaWVudGF0aW9uL25vZGVfbW9kdWxlcy9yb2J1c3Qtc3VidHJhY3Qvcm9idXN0LWRpZmYuanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL2FudHMvQW50Q29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvaW5jcmVtZW50YWwtY29udmV4LWh1bGwvbm9kZV9tb2R1bGVzL3JvYnVzdC1vcmllbnRhdGlvbi9ub2RlX21vZHVsZXMvcm9idXN0LXN1bS9yb2J1c3Qtc3VtLmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9hbnRzL0FudENvbG9ueS9ub2RlX21vZHVsZXMvZGVsYXVuYXktdHJpYW5ndWxhdGUvbm9kZV9tb2R1bGVzL2luY3JlbWVudGFsLWNvbnZleC1odWxsL25vZGVfbW9kdWxlcy9yb2J1c3Qtb3JpZW50YXRpb24vbm9kZV9tb2R1bGVzL3R3by1wcm9kdWN0L3R3by1wcm9kdWN0LmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9hbnRzL0FudENvbG9ueS9ub2RlX21vZHVsZXMvZGVsYXVuYXktdHJpYW5ndWxhdGUvbm9kZV9tb2R1bGVzL2luY3JlbWVudGFsLWNvbnZleC1odWxsL25vZGVfbW9kdWxlcy9yb2J1c3Qtb3JpZW50YXRpb24vb3JpZW50YXRpb24uanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL2FudHMvQW50Q29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvaW5jcmVtZW50YWwtY29udmV4LWh1bGwvbm9kZV9tb2R1bGVzL3NpbXBsaWNpYWwtY29tcGxleC9ub2RlX21vZHVsZXMvYml0LXR3aWRkbGUvdHdpZGRsZS5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vYW50cy9BbnRDb2xvbnkvbm9kZV9tb2R1bGVzL2RlbGF1bmF5LXRyaWFuZ3VsYXRlL25vZGVfbW9kdWxlcy9pbmNyZW1lbnRhbC1jb252ZXgtaHVsbC9ub2RlX21vZHVsZXMvc2ltcGxpY2lhbC1jb21wbGV4L25vZGVfbW9kdWxlcy91bmlvbi1maW5kL2luZGV4LmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9hbnRzL0FudENvbG9ueS9ub2RlX21vZHVsZXMvZGVsYXVuYXktdHJpYW5ndWxhdGUvbm9kZV9tb2R1bGVzL2luY3JlbWVudGFsLWNvbnZleC1odWxsL25vZGVfbW9kdWxlcy9zaW1wbGljaWFsLWNvbXBsZXgvdG9wb2xvZ3kuanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL2FudHMvQW50Q29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvdW5pcS91bmlxLmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9hbnRzL0FudENvbG9ueS9ub2RlX21vZHVsZXMvZGVsYXVuYXktdHJpYW5ndWxhdGUvdHJpYW5ndWxhdGUuanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL2FudHMvQW50Q29sb255L25vZGVfbW9kdWxlcy9wYXJzZS1zdmctcGF0aC9pbmRleC5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vYW50cy9BbnRDb2xvbnkvc3JjL2FudC5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vYW50cy9BbnRDb2xvbnkvc3JjL2FudHNHcm91cC5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vYW50cy9BbnRDb2xvbnkvc3JjL2NyZWF0ZUVkZ2VzLmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9hbnRzL0FudENvbG9ueS9zcmMvZWRnZS5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vYW50cy9BbnRDb2xvbnkvc3JjL2luaXRpYWxpemVQb2ludHMuanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL2FudHMvQW50Q29sb255L3NyYy9tb3VzZS5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vYW50cy9BbnRDb2xvbnkvc3JjL3BvaW50LmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9hbnRzL0FudENvbG9ueS9zcmMvcmVuZGVyaW5nLmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9hbnRzL0FudENvbG9ueS9zcmMvc3ZnUGF0aC5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vYW50cy9BbnRDb2xvbnkvc3JjL3V0aWxpdGllcy5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vYW50cy9BbnRDb2xvbnkvc3JjL3ZlY3Rvci5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vYW50cy9BbnRDb2xvbnkvc3RhcnQuanMiXSwibmFtZXMiOltdLCJtYXBwaW5ncyI6IkFBQUE7QUNBQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUN2Q0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUM3YkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNoQkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNqREE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQzNKQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDM0pBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNoQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDN0xBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQzVNQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUN6REE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDdFZBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ3pEQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDOUpBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDdkRBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDbk5BO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNuQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNyRkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDOUVBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDNUtBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDbkJBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ1RBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQzdPQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ2hCQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQzVCQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ25CQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSIsImZpbGUiOiJnZW5lcmF0ZWQuanMiLCJzb3VyY2VSb290IjoiIiwic291cmNlc0NvbnRlbnQiOlsiKGZ1bmN0aW9uIGUodCxuLHIpe2Z1bmN0aW9uIHMobyx1KXtpZighbltvXSl7aWYoIXRbb10pe3ZhciBhPXR5cGVvZiByZXF1aXJlPT1cImZ1bmN0aW9uXCImJnJlcXVpcmU7aWYoIXUmJmEpcmV0dXJuIGEobywhMCk7aWYoaSlyZXR1cm4gaShvLCEwKTt2YXIgZj1uZXcgRXJyb3IoXCJDYW5ub3QgZmluZCBtb2R1bGUgJ1wiK28rXCInXCIpO3Rocm93IGYuY29kZT1cIk1PRFVMRV9OT1RfRk9VTkRcIixmfXZhciBsPW5bb109e2V4cG9ydHM6e319O3Rbb11bMF0uY2FsbChsLmV4cG9ydHMsZnVuY3Rpb24oZSl7dmFyIG49dFtvXVsxXVtlXTtyZXR1cm4gcyhuP246ZSl9LGwsbC5leHBvcnRzLGUsdCxuLHIpfXJldHVybiBuW29dLmV4cG9ydHN9dmFyIGk9dHlwZW9mIHJlcXVpcmU9PVwiZnVuY3Rpb25cIiYmcmVxdWlyZTtmb3IodmFyIG89MDtvPHIubGVuZ3RoO28rKylzKHJbb10pO3JldHVybiBzfSkiLCIndXNlIHN0cmljdCc7XG5cbnZhciBpbml0UmVuZGVyaW5nID0gcmVxdWlyZSgnLi9zcmMvcmVuZGVyaW5nLmpzJyk7XG52YXIgaW5pdGlhbGl6ZVBvaW50cyA9IHJlcXVpcmUoJy4vc3JjL2luaXRpYWxpemVQb2ludHMuanMnKTtcbnZhciBjcmVhdGVFZGdlcyA9IHJlcXVpcmUoJy4vc3JjL2NyZWF0ZUVkZ2VzLmpzJyk7XG4vLyB2YXIgaW5pdEFudHMgPSByZXF1aXJlKCcuL3NyYy9pbml0aWFsaXplQW50cycpO1xuXG5tb2R1bGUuZXhwb3J0cyA9IGZ1bmN0aW9uIGluaXQoY29udGFpbmVyRWxlbWVudCwgb3B0aW9ucyl7XG5cblx0dmFyIHJlbmRlciwgcG9pbnRzSW5mb3MsIGVkZ2VzLCBwb3B1bGF0aW9uLCBwb2ludHNNYXA7XG5cblxuXHRmdW5jdGlvbiBfaW5pdChjb250YWluZXJFbGVtZW50LCBvcHRpb25zKXtcblx0XHRwb2ludHNJbmZvcyA9IGluaXRpYWxpemVQb2ludHMob3B0aW9ucy5uYlN0YXJ0LCBvcHRpb25zLm5iUmFuZCk7XG5cdFx0ZWRnZXMgPSBjcmVhdGVFZGdlcyhwb2ludHNJbmZvcy5wb2ludHMpO1xuXHRcdC8vIHBvcHVsYXRpb24gPSBvcHRpb25zLm5iQW50cztcblx0XHQvLyBwb3B1bGF0aW9uID0gaW5pdEFudHMoY29udGFpbmVyRWxlbWVudCwgcG9pbnRzSW5mb3MsIG9wdGlvbnMpO1xuXHRcdHBvaW50c01hcCA9IHtcblx0XHRcdHBvaW50c0luZm9zOiBwb2ludHNJbmZvcyxcblx0XHRcdGVkZ2VzOiBlZGdlc1xuXHRcdFx0Ly8gcG9wdWxhdGlvbjogcG9wdWxhdGlvblxuXHRcdH07XG5cdFx0cmVuZGVyID0gaW5pdFJlbmRlcmluZyhjb250YWluZXJFbGVtZW50LCBwb2ludHNNYXAsIG9wdGlvbnMpO1xuXHR9XG5cblx0X2luaXQoY29udGFpbmVyRWxlbWVudCwgb3B0aW9ucyk7XG5cblx0cmV0dXJuIHtcblx0XHR0b2dnbGVQbGF5UGF1c2U6IGZ1bmN0aW9uKCl7IHJlbmRlci50b2dnbGVQbGF5UGF1c2UoKSB9LFxuXHRcdGNoYW5nZU9wdGlvbnM6IGZ1bmN0aW9uKG9wdHMpe1xuXHRcdFx0cmVuZGVyLm1vZGlmeUFudHMob3B0cyk7XG5cdFx0fSxcblx0XHRyZXNldDogZnVuY3Rpb24ob3B0cyl7XG5cdFx0XHRyZW5kZXIucmVzZXQoKTtcblxuXHRcdFx0XHQvLyByZXNldCBlbGVtZW50c1xuXHRcdFx0X2luaXQoY29udGFpbmVyRWxlbWVudCwgb3B0cyk7XG5cdFx0fVxuXHR9O1xufTsiLCJcInVzZSBzdHJpY3RcIlxuXG4vL0hpZ2ggbGV2ZWwgaWRlYTpcbi8vIDEuIFVzZSBDbGFya3NvbidzIGluY3JlbWVudGFsIGNvbnN0cnVjdGlvbiB0byBmaW5kIGNvbnZleCBodWxsXG4vLyAyLiBQb2ludCBsb2NhdGlvbiBpbiB0cmlhbmd1bGF0aW9uIGJ5IGp1bXAgYW5kIHdhbGtcblxubW9kdWxlLmV4cG9ydHMgPSBpbmNyZW1lbnRhbENvbnZleEh1bGxcblxudmFyIG9yaWVudCA9IHJlcXVpcmUoXCJyb2J1c3Qtb3JpZW50YXRpb25cIilcbnZhciBjb21wYXJlQ2VsbCA9IHJlcXVpcmUoXCJzaW1wbGljaWFsLWNvbXBsZXhcIikuY29tcGFyZUNlbGxzXG5cbmZ1bmN0aW9uIGNvbXBhcmVJbnQoYSwgYikge1xuICByZXR1cm4gYSAtIGJcbn1cblxuZnVuY3Rpb24gU2ltcGxleCh2ZXJ0aWNlcywgYWRqYWNlbnQsIGJvdW5kYXJ5KSB7XG4gIHRoaXMudmVydGljZXMgPSB2ZXJ0aWNlc1xuICB0aGlzLmFkamFjZW50ID0gYWRqYWNlbnRcbiAgdGhpcy5ib3VuZGFyeSA9IGJvdW5kYXJ5XG4gIHRoaXMubGFzdFZpc2l0ZWQgPSAtMVxufVxuXG5TaW1wbGV4LnByb3RvdHlwZS5mbGlwID0gZnVuY3Rpb24oKSB7XG4gIHZhciB0ID0gdGhpcy52ZXJ0aWNlc1swXVxuICB0aGlzLnZlcnRpY2VzWzBdID0gdGhpcy52ZXJ0aWNlc1sxXVxuICB0aGlzLnZlcnRpY2VzWzFdID0gdFxuICB2YXIgdSA9IHRoaXMuYWRqYWNlbnRbMF1cbiAgdGhpcy5hZGphY2VudFswXSA9IHRoaXMuYWRqYWNlbnRbMV1cbiAgdGhpcy5hZGphY2VudFsxXSA9IHVcbn1cblxuZnVuY3Rpb24gR2x1ZUZhY2V0KHZlcnRpY2VzLCBjZWxsLCBpbmRleCkge1xuICB0aGlzLnZlcnRpY2VzID0gdmVydGljZXNcbiAgdGhpcy5jZWxsID0gY2VsbFxuICB0aGlzLmluZGV4ID0gaW5kZXhcbn1cblxuZnVuY3Rpb24gY29tcGFyZUdsdWUoYSwgYikge1xuICByZXR1cm4gY29tcGFyZUNlbGwoYS52ZXJ0aWNlcywgYi52ZXJ0aWNlcylcbn1cblxuZnVuY3Rpb24gYmFrZU9yaWVudChkKSB7XG4gIHZhciBjb2RlID0gW1wiZnVuY3Rpb24gb3JpZW50KCl7dmFyIHR1cGxlPXRoaXMudHVwbGU7cmV0dXJuIHRlc3QoXCJdXG4gIGZvcih2YXIgaT0wOyBpPD1kOyArK2kpIHtcbiAgICBpZihpID4gMCkge1xuICAgICAgY29kZS5wdXNoKFwiLFwiKVxuICAgIH1cbiAgICBjb2RlLnB1c2goXCJ0dXBsZVtcIiwgaSwgXCJdXCIpXG4gIH1cbiAgY29kZS5wdXNoKFwiKX1yZXR1cm4gb3JpZW50XCIpXG4gIHZhciBwcm9jID0gbmV3IEZ1bmN0aW9uKFwidGVzdFwiLCBjb2RlLmpvaW4oXCJcIikpXG4gIHZhciB0ZXN0ID0gb3JpZW50W2QrMV1cbiAgaWYoIXRlc3QpIHtcbiAgICB0ZXN0ID0gb3JpZW50XG4gIH1cbiAgcmV0dXJuIHByb2ModGVzdClcbn1cblxudmFyIEJBS0VEID0gW11cblxuZnVuY3Rpb24gVHJpYW5ndWxhdGlvbihkaW1lbnNpb24sIHZlcnRpY2VzLCBzaW1wbGljZXMpIHtcbiAgdGhpcy5kaW1lbnNpb24gPSBkaW1lbnNpb25cbiAgdGhpcy52ZXJ0aWNlcyA9IHZlcnRpY2VzXG4gIHRoaXMuc2ltcGxpY2VzID0gc2ltcGxpY2VzXG4gIHRoaXMuaW50ZXJpb3IgPSBzaW1wbGljZXMuZmlsdGVyKGZ1bmN0aW9uKGMpIHtcbiAgICByZXR1cm4gIWMuYm91bmRhcnlcbiAgfSlcblxuICB0aGlzLnR1cGxlID0gbmV3IEFycmF5KGRpbWVuc2lvbisxKVxuICBmb3IodmFyIGk9MDsgaTw9ZGltZW5zaW9uOyArK2kpIHtcbiAgICB0aGlzLnR1cGxlW2ldID0gdGhpcy52ZXJ0aWNlc1tpXVxuICB9XG5cbiAgdmFyIG8gPSBCQUtFRFtkaW1lbnNpb25dXG4gIGlmKCFvKSB7XG4gICAgbyA9IEJBS0VEW2RpbWVuc2lvbl0gPSBiYWtlT3JpZW50KGRpbWVuc2lvbilcbiAgfVxuICB0aGlzLm9yaWVudCA9IG9cbn1cblxudmFyIHByb3RvID0gVHJpYW5ndWxhdGlvbi5wcm90b3R5cGVcblxuLy9EZWdlbmVyYXRlIHNpdHVhdGlvbiB3aGVyZSB3ZSBhcmUgb24gYm91bmRhcnksIGJ1dCBjb3BsYW5hciB0byBmYWNlXG5wcm90by5oYW5kbGVCb3VuZGFyeURlZ2VuZXJhY3kgPSBmdW5jdGlvbihjZWxsLCBwb2ludCkge1xuICB2YXIgZCA9IHRoaXMuZGltZW5zaW9uXG4gIHZhciBuID0gdGhpcy52ZXJ0aWNlcy5sZW5ndGggLSAxXG4gIHZhciB0dXBsZSA9IHRoaXMudHVwbGVcbiAgdmFyIHZlcnRzID0gdGhpcy52ZXJ0aWNlc1xuXG4gIC8vRHVtYiBzb2x1dGlvbjogSnVzdCBkbyBkZnMgZnJvbSBib3VuZGFyeSBjZWxsIHVudGlsIHdlIGZpbmQgYW55IHBlYWssIG9yIHRlcm1pbmF0ZVxuICB2YXIgdG9WaXNpdCA9IFsgY2VsbCBdXG4gIGNlbGwubGFzdFZpc2l0ZWQgPSAtblxuICB3aGlsZSh0b1Zpc2l0Lmxlbmd0aCA+IDApIHtcbiAgICBjZWxsID0gdG9WaXNpdC5wb3AoKVxuICAgIHZhciBjZWxsVmVydHMgPSBjZWxsLnZlcnRpY2VzXG4gICAgdmFyIGNlbGxBZGogPSBjZWxsLmFkamFjZW50XG4gICAgZm9yKHZhciBpPTA7IGk8PWQ7ICsraSkge1xuICAgICAgdmFyIG5laWdoYm9yID0gY2VsbEFkaltpXVxuICAgICAgaWYoIW5laWdoYm9yLmJvdW5kYXJ5IHx8IG5laWdoYm9yLmxhc3RWaXNpdGVkIDw9IC1uKSB7XG4gICAgICAgIGNvbnRpbnVlXG4gICAgICB9XG4gICAgICB2YXIgbnYgPSBuZWlnaGJvci52ZXJ0aWNlc1xuICAgICAgZm9yKHZhciBqPTA7IGo8PWQ7ICsraikge1xuICAgICAgICB2YXIgdnYgPSBudltqXVxuICAgICAgICBpZih2diA8IDApIHtcbiAgICAgICAgICB0dXBsZVtqXSA9IHBvaW50XG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgdHVwbGVbal0gPSB2ZXJ0c1t2dl1cbiAgICAgICAgfVxuICAgICAgfVxuICAgICAgdmFyIG8gPSB0aGlzLm9yaWVudCgpXG4gICAgICBpZihvID4gMCkge1xuICAgICAgICByZXR1cm4gbmVpZ2hib3JcbiAgICAgIH1cbiAgICAgIG5laWdoYm9yLmxhc3RWaXNpdGVkID0gLW5cbiAgICAgIGlmKG8gPT09IDApIHtcbiAgICAgICAgdG9WaXNpdC5wdXNoKG5laWdoYm9yKVxuICAgICAgfVxuICAgIH1cbiAgfVxuICByZXR1cm4gbnVsbFxufVxuXG5wcm90by53YWxrID0gZnVuY3Rpb24ocG9pbnQsIHJhbmRvbSkge1xuICAvL0FsaWFzIGxvY2FsIHByb3BlcnRpZXNcbiAgdmFyIG4gPSB0aGlzLnZlcnRpY2VzLmxlbmd0aCAtIDFcbiAgdmFyIGQgPSB0aGlzLmRpbWVuc2lvblxuICB2YXIgdmVydHMgPSB0aGlzLnZlcnRpY2VzXG4gIHZhciB0dXBsZSA9IHRoaXMudHVwbGVcblxuICAvL0NvbXB1dGUgaW5pdGlhbCBqdW1wIGNlbGxcbiAgdmFyIGluaXRJbmRleCA9IHJhbmRvbSA/ICh0aGlzLmludGVyaW9yLmxlbmd0aCAqIE1hdGgucmFuZG9tKCkpfDAgOiAodGhpcy5pbnRlcmlvci5sZW5ndGgtMSlcbiAgdmFyIGNlbGwgPSB0aGlzLmludGVyaW9yWyBpbml0SW5kZXggXVxuXG4gIC8vU3RhcnQgd2Fsa2luZ1xub3V0ZXJMb29wOlxuICB3aGlsZSghY2VsbC5ib3VuZGFyeSkge1xuICAgIHZhciBjZWxsVmVydHMgPSBjZWxsLnZlcnRpY2VzXG4gICAgdmFyIGNlbGxBZGogPSBjZWxsLmFkamFjZW50XG5cbiAgICBmb3IodmFyIGk9MDsgaTw9ZDsgKytpKSB7XG4gICAgICB0dXBsZVtpXSA9IHZlcnRzW2NlbGxWZXJ0c1tpXV1cbiAgICB9XG4gICAgY2VsbC5sYXN0VmlzaXRlZCA9IG5cblxuICAgIC8vRmluZCBmYXJ0aGVzdCBhZGphY2VudCBjZWxsXG4gICAgZm9yKHZhciBpPTA7IGk8PWQ7ICsraSkge1xuICAgICAgdmFyIG5laWdoYm9yID0gY2VsbEFkaltpXVxuICAgICAgaWYobmVpZ2hib3IubGFzdFZpc2l0ZWQgPj0gbikge1xuICAgICAgICBjb250aW51ZVxuICAgICAgfVxuICAgICAgdmFyIHByZXYgPSB0dXBsZVtpXVxuICAgICAgdHVwbGVbaV0gPSBwb2ludFxuICAgICAgdmFyIG8gPSB0aGlzLm9yaWVudCgpXG4gICAgICB0dXBsZVtpXSA9IHByZXZcbiAgICAgIGlmKG8gPCAwKSB7XG4gICAgICAgIGNlbGwgPSBuZWlnaGJvclxuICAgICAgICBjb250aW51ZSBvdXRlckxvb3BcbiAgICAgIH0gZWxzZSB7XG4gICAgICAgIGlmKCFuZWlnaGJvci5ib3VuZGFyeSkge1xuICAgICAgICAgIG5laWdoYm9yLmxhc3RWaXNpdGVkID0gblxuICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgIG5laWdoYm9yLmxhc3RWaXNpdGVkID0gLW5cbiAgICAgICAgfVxuICAgICAgfVxuICAgIH1cbiAgICByZXR1cm5cbiAgfVxuXG4gIHJldHVybiBjZWxsXG59XG5cbnByb3RvLmFkZFBlYWtzID0gZnVuY3Rpb24ocG9pbnQsIGNlbGwpIHtcbiAgdmFyIG4gPSB0aGlzLnZlcnRpY2VzLmxlbmd0aCAtIDFcbiAgdmFyIGQgPSB0aGlzLmRpbWVuc2lvblxuICB2YXIgdmVydHMgPSB0aGlzLnZlcnRpY2VzXG4gIHZhciB0dXBsZSA9IHRoaXMudHVwbGVcbiAgdmFyIGludGVyaW9yID0gdGhpcy5pbnRlcmlvclxuICB2YXIgc2ltcGxpY2VzID0gdGhpcy5zaW1wbGljZXNcblxuICAvL1dhbGtpbmcgZmluaXNoZWQgYXQgYm91bmRhcnksIHRpbWUgdG8gYWRkIHBlYWtzXG4gIHZhciB0b3Zpc2l0ID0gWyBjZWxsIF1cblxuICAvL1N0cmV0Y2ggaW5pdGlhbCBib3VuZGFyeSBjZWxsIGludG8gYSBwZWFrXG4gIGNlbGwubGFzdFZpc2l0ZWQgPSBuXG4gIGNlbGwudmVydGljZXNbY2VsbC52ZXJ0aWNlcy5pbmRleE9mKC0xKV0gPSBuXG4gIGNlbGwuYm91bmRhcnkgPSBmYWxzZVxuICBpbnRlcmlvci5wdXNoKGNlbGwpXG5cbiAgLy9SZWNvcmQgYSBsaXN0IG9mIGFsbCBuZXcgYm91bmRhcmllcyBjcmVhdGVkIGJ5IGFkZGVkIHBlYWtzIHNvIHdlIGNhbiBnbHVlIHRoZW0gdG9nZXRoZXIgd2hlbiB3ZSBhcmUgYWxsIGRvbmVcbiAgdmFyIGdsdWVGYWNldHMgPSBbXVxuXG4gIC8vRG8gYSB0cmF2ZXJzYWwgb2YgdGhlIGJvdW5kYXJ5IHdhbGtpbmcgb3V0d2FyZCBmcm9tIHN0YXJ0aW5nIHBlYWtcbiAgd2hpbGUodG92aXNpdC5sZW5ndGggPiAwKSB7XG4gICAgLy9Qb3Agb2ZmIHBlYWsgYW5kIHdhbGsgb3ZlciBhZGphY2VudCBjZWxsc1xuICAgIHZhciBjZWxsID0gdG92aXNpdC5wb3AoKVxuICAgIHZhciBjZWxsVmVydHMgPSBjZWxsLnZlcnRpY2VzXG4gICAgdmFyIGNlbGxBZGogPSBjZWxsLmFkamFjZW50XG4gICAgdmFyIGluZGV4T2ZOID0gY2VsbFZlcnRzLmluZGV4T2YobilcbiAgICBpZihpbmRleE9mTiA8IDApIHtcbiAgICAgIGNvbnRpbnVlXG4gICAgfVxuXG4gICAgZm9yKHZhciBpPTA7IGk8PWQ7ICsraSkge1xuICAgICAgaWYoaSA9PT0gaW5kZXhPZk4pIHtcbiAgICAgICAgY29udGludWVcbiAgICAgIH1cblxuICAgICAgLy9Gb3IgZWFjaCBib3VuZGFyeSBuZWlnaGJvciBvZiB0aGUgY2VsbFxuICAgICAgdmFyIG5laWdoYm9yID0gY2VsbEFkaltpXVxuICAgICAgaWYoIW5laWdoYm9yLmJvdW5kYXJ5IHx8IG5laWdoYm9yLmxhc3RWaXNpdGVkID49IG4pIHtcbiAgICAgICAgY29udGludWVcbiAgICAgIH1cblxuICAgICAgdmFyIG52ID0gbmVpZ2hib3IudmVydGljZXNcblxuICAgICAgLy9UZXN0IGlmIG5laWdoYm9yIGlzIGEgcGVha1xuICAgICAgaWYobmVpZ2hib3IubGFzdFZpc2l0ZWQgIT09IC1uKSB7ICAgICAgXG4gICAgICAgIC8vQ29tcHV0ZSBvcmllbnRhdGlvbiBvZiBwIHJlbGF0aXZlIHRvIGVhY2ggYm91bmRhcnkgcGVha1xuICAgICAgICB2YXIgaW5kZXhPZk5lZzEgPSAwXG4gICAgICAgIGZvcih2YXIgaj0wOyBqPD1kOyArK2opIHtcbiAgICAgICAgICBpZihudltqXSA8IDApIHtcbiAgICAgICAgICAgIGluZGV4T2ZOZWcxID0galxuICAgICAgICAgICAgdHVwbGVbal0gPSBwb2ludFxuICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICB0dXBsZVtqXSA9IHZlcnRzW252W2pdXVxuICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgICAgICB2YXIgbyA9IHRoaXMub3JpZW50KClcblxuICAgICAgICAvL1Rlc3QgaWYgbmVpZ2hib3IgY2VsbCBpcyBhbHNvIGEgcGVha1xuICAgICAgICBpZihvID4gMCkge1xuICAgICAgICAgIG52W2luZGV4T2ZOZWcxXSA9IG5cbiAgICAgICAgICBuZWlnaGJvci5ib3VuZGFyeSA9IGZhbHNlXG4gICAgICAgICAgaW50ZXJpb3IucHVzaChuZWlnaGJvcilcbiAgICAgICAgICB0b3Zpc2l0LnB1c2gobmVpZ2hib3IpXG4gICAgICAgICAgbmVpZ2hib3IubGFzdFZpc2l0ZWQgPSBuXG4gICAgICAgICAgY29udGludWVcbiAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICBuZWlnaGJvci5sYXN0VmlzaXRlZCA9IC1uXG4gICAgICAgIH1cbiAgICAgIH1cblxuICAgICAgdmFyIG5hID0gbmVpZ2hib3IuYWRqYWNlbnRcblxuICAgICAgLy9PdGhlcndpc2UsIHJlcGxhY2UgbmVpZ2hib3Igd2l0aCBuZXcgZmFjZVxuICAgICAgdmFyIHZ2ZXJ0cyA9IGNlbGxWZXJ0cy5zbGljZSgpXG4gICAgICB2YXIgdmFkaiA9IGNlbGxBZGouc2xpY2UoKVxuICAgICAgdmFyIG5jZWxsID0gbmV3IFNpbXBsZXgodnZlcnRzLCB2YWRqLCB0cnVlKVxuICAgICAgc2ltcGxpY2VzLnB1c2gobmNlbGwpXG5cbiAgICAgIC8vQ29ubmVjdCB0byBuZWlnaGJvclxuICAgICAgdmFyIG9wcG9zaXRlID0gbmEuaW5kZXhPZihjZWxsKVxuICAgICAgaWYob3Bwb3NpdGUgPCAwKSB7XG4gICAgICAgIGNvbnRpbnVlXG4gICAgICB9XG4gICAgICBuYVtvcHBvc2l0ZV0gPSBuY2VsbFxuICAgICAgdmFkaltpbmRleE9mTl0gPSBuZWlnaGJvclxuXG4gICAgICAvL0Nvbm5lY3QgdG8gY2VsbFxuICAgICAgdnZlcnRzW2ldID0gLTFcbiAgICAgIHZhZGpbaV0gPSBjZWxsXG4gICAgICBjZWxsQWRqW2ldID0gbmNlbGxcblxuICAgICAgLy9GbGlwIGZhY2V0XG4gICAgICBuY2VsbC5mbGlwKClcblxuICAgICAgLy9BZGQgdG8gZ2x1ZSBsaXN0XG4gICAgICBmb3IodmFyIGo9MDsgajw9ZDsgKytqKSB7XG4gICAgICAgIHZhciB1dSA9IHZ2ZXJ0c1tqXVxuICAgICAgICBpZih1dSA8IDAgfHwgdXUgPT09IG4pIHtcbiAgICAgICAgICBjb250aW51ZVxuICAgICAgICB9XG4gICAgICAgIHZhciBuZmFjZSA9IG5ldyBBcnJheShkLTEpXG4gICAgICAgIHZhciBucHRyID0gMFxuICAgICAgICBmb3IodmFyIGs9MDsgazw9ZDsgKytrKSB7XG4gICAgICAgICAgdmFyIHZ2ID0gdnZlcnRzW2tdXG4gICAgICAgICAgaWYodnYgPCAwIHx8IGsgPT09IGopIHtcbiAgICAgICAgICAgIGNvbnRpbnVlXG4gICAgICAgICAgfVxuICAgICAgICAgIG5mYWNlW25wdHIrK10gPSB2dlxuICAgICAgICB9XG4gICAgICAgIGdsdWVGYWNldHMucHVzaChuZXcgR2x1ZUZhY2V0KG5mYWNlLCBuY2VsbCwgaikpXG4gICAgICB9XG4gICAgfVxuICB9XG5cbiAgLy9HbHVlIGJvdW5kYXJ5IGZhY2V0cyB0b2dldGhlclxuICBnbHVlRmFjZXRzLnNvcnQoY29tcGFyZUdsdWUpXG5cbiAgZm9yKHZhciBpPTA7IGkrMTxnbHVlRmFjZXRzLmxlbmd0aDsgaSs9Mikge1xuICAgIHZhciBhID0gZ2x1ZUZhY2V0c1tpXVxuICAgIHZhciBiID0gZ2x1ZUZhY2V0c1tpKzFdXG4gICAgdmFyIGFpID0gYS5pbmRleFxuICAgIHZhciBiaSA9IGIuaW5kZXhcbiAgICBpZihhaSA8IDAgfHwgYmkgPCAwKSB7XG4gICAgICBjb250aW51ZVxuICAgIH1cbiAgICBhLmNlbGwuYWRqYWNlbnRbYS5pbmRleF0gPSBiLmNlbGxcbiAgICBiLmNlbGwuYWRqYWNlbnRbYi5pbmRleF0gPSBhLmNlbGxcbiAgfVxufVxuXG5wcm90by5pbnNlcnQgPSBmdW5jdGlvbihwb2ludCwgcmFuZG9tKSB7XG4gIC8vQWRkIHBvaW50XG4gIHZhciB2ZXJ0cyA9IHRoaXMudmVydGljZXNcbiAgdmVydHMucHVzaChwb2ludClcblxuICB2YXIgY2VsbCA9IHRoaXMud2Fsayhwb2ludCwgcmFuZG9tKVxuICBpZighY2VsbCkge1xuICAgIHJldHVyblxuICB9XG5cbiAgLy9BbGlhcyBsb2NhbCBwcm9wZXJ0aWVzXG4gIHZhciBkID0gdGhpcy5kaW1lbnNpb25cbiAgdmFyIHR1cGxlID0gdGhpcy50dXBsZVxuXG4gIC8vRGVnZW5lcmF0ZSBjYXNlOiBJZiBwb2ludCBpcyBjb3BsYW5hciB0byBjZWxsLCB0aGVuIHdhbGsgdW50aWwgd2UgZmluZCBhIG5vbi1kZWdlbmVyYXRlIGJvdW5kYXJ5XG4gIGZvcih2YXIgaT0wOyBpPD1kOyArK2kpIHtcbiAgICB2YXIgdnYgPSBjZWxsLnZlcnRpY2VzW2ldXG4gICAgaWYodnYgPCAwKSB7XG4gICAgICB0dXBsZVtpXSA9IHBvaW50XG4gICAgfSBlbHNlIHtcbiAgICAgIHR1cGxlW2ldID0gdmVydHNbdnZdXG4gICAgfVxuICB9XG4gIHZhciBvID0gdGhpcy5vcmllbnQodHVwbGUpXG4gIGlmKG8gPCAwKSB7XG4gICAgcmV0dXJuXG4gIH0gZWxzZSBpZihvID09PSAwKSB7XG4gICAgY2VsbCA9IHRoaXMuaGFuZGxlQm91bmRhcnlEZWdlbmVyYWN5KGNlbGwsIHBvaW50KVxuICAgIGlmKCFjZWxsKSB7XG4gICAgICByZXR1cm5cbiAgICB9XG4gIH1cblxuICAvL0FkZCBwZWFrc1xuICB0aGlzLmFkZFBlYWtzKHBvaW50LCBjZWxsKVxufVxuXG4vL0V4dHJhY3QgYWxsIGJvdW5kYXJ5IGNlbGxzXG5wcm90by5ib3VuZGFyeSA9IGZ1bmN0aW9uKCkge1xuICB2YXIgZCA9IHRoaXMuZGltZW5zaW9uXG4gIHZhciBib3VuZGFyeSA9IFtdXG4gIHZhciBjZWxscyA9IHRoaXMuc2ltcGxpY2VzXG4gIHZhciBuYyA9IGNlbGxzLmxlbmd0aFxuICBmb3IodmFyIGk9MDsgaTxuYzsgKytpKSB7XG4gICAgdmFyIGMgPSBjZWxsc1tpXVxuICAgIGlmKGMuYm91bmRhcnkpIHtcbiAgICAgIHZhciBiY2VsbCA9IG5ldyBBcnJheShkKVxuICAgICAgdmFyIGN2ID0gYy52ZXJ0aWNlc1xuICAgICAgdmFyIHB0ciA9IDBcbiAgICAgIHZhciBwYXJpdHkgPSAwXG4gICAgICBmb3IodmFyIGo9MDsgajw9ZDsgKytqKSB7XG4gICAgICAgIGlmKGN2W2pdID49IDApIHtcbiAgICAgICAgICBiY2VsbFtwdHIrK10gPSBjdltqXVxuICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgIHBhcml0eSA9IGomMVxuICAgICAgICB9XG4gICAgICB9XG4gICAgICBpZihwYXJpdHkgPT09IChkJjEpKSB7XG4gICAgICAgIHZhciB0ID0gYmNlbGxbMF1cbiAgICAgICAgYmNlbGxbMF0gPSBiY2VsbFsxXVxuICAgICAgICBiY2VsbFsxXSA9IHRcbiAgICAgIH1cbiAgICAgIGJvdW5kYXJ5LnB1c2goYmNlbGwpXG4gICAgfVxuICB9XG4gIHJldHVybiBib3VuZGFyeVxufVxuXG5mdW5jdGlvbiBpbmNyZW1lbnRhbENvbnZleEh1bGwocG9pbnRzLCByYW5kb21TZWFyY2gpIHtcbiAgdmFyIG4gPSBwb2ludHMubGVuZ3RoXG4gIGlmKG4gPT09IDApIHtcbiAgICB0aHJvdyBuZXcgRXJyb3IoXCJNdXN0IGhhdmUgYXQgbGVhc3QgZCsxIHBvaW50c1wiKVxuICB9XG4gIHZhciBkID0gcG9pbnRzWzBdLmxlbmd0aFxuICBpZihuIDw9IGQpIHtcbiAgICB0aHJvdyBuZXcgRXJyb3IoXCJNdXN0IGlucHV0IGF0IGxlYXN0IGQrMSBwb2ludHNcIilcbiAgfVxuXG4gIC8vRklYTUU6IFRoaXMgY291bGQgYmUgZGVnZW5lcmF0ZSwgYnV0IG5lZWQgdG8gc2VsZWN0IGQrMSBub24tY29wbGFuYXIgcG9pbnRzIHRvIGJvb3RzdHJhcCBwcm9jZXNzXG4gIHZhciBpbml0aWFsU2ltcGxleCA9IHBvaW50cy5zbGljZSgwLCBkKzEpXG5cbiAgLy9NYWtlIHN1cmUgaW5pdGlhbCBzaW1wbGV4IGlzIHBvc2l0aXZlbHkgb3JpZW50ZWRcbiAgdmFyIG8gPSBvcmllbnQuYXBwbHkodm9pZCAwLCBpbml0aWFsU2ltcGxleClcbiAgaWYobyA9PT0gMCkge1xuICAgIHRocm93IG5ldyBFcnJvcihcIklucHV0IG5vdCBpbiBnZW5lcmFsIHBvc2l0aW9uXCIpXG4gIH1cbiAgdmFyIGluaXRpYWxDb29yZHMgPSBuZXcgQXJyYXkoZCsxKVxuICBmb3IodmFyIGk9MDsgaTw9ZDsgKytpKSB7XG4gICAgaW5pdGlhbENvb3Jkc1tpXSA9IGlcbiAgfVxuICBpZihvIDwgMCkge1xuICAgIGluaXRpYWxDb29yZHNbMF0gPSAxXG4gICAgaW5pdGlhbENvb3Jkc1sxXSA9IDBcbiAgfVxuXG4gIC8vQ3JlYXRlIGluaXRpYWwgdG9wb2xvZ2ljYWwgaW5kZXgsIGdsdWUgcG9pbnRlcnMgdG9nZXRoZXIgKGtpbmQgb2YgbWVzc3kpXG4gIHZhciBpbml0aWFsQ2VsbCA9IG5ldyBTaW1wbGV4KGluaXRpYWxDb29yZHMsIG5ldyBBcnJheShkKzEpLCBmYWxzZSlcbiAgdmFyIGJvdW5kYXJ5ID0gaW5pdGlhbENlbGwuYWRqYWNlbnRcbiAgdmFyIGxpc3QgPSBuZXcgQXJyYXkoZCsyKVxuICBmb3IodmFyIGk9MDsgaTw9ZDsgKytpKSB7XG4gICAgdmFyIHZlcnRzID0gaW5pdGlhbENvb3Jkcy5zbGljZSgpXG4gICAgZm9yKHZhciBqPTA7IGo8PWQ7ICsraikge1xuICAgICAgaWYoaiA9PT0gaSkge1xuICAgICAgICB2ZXJ0c1tqXSA9IC0xXG4gICAgICB9XG4gICAgfVxuICAgIHZhciB0ID0gdmVydHNbMF1cbiAgICB2ZXJ0c1swXSA9IHZlcnRzWzFdXG4gICAgdmVydHNbMV0gPSB0XG4gICAgdmFyIGNlbGwgPSBuZXcgU2ltcGxleCh2ZXJ0cywgbmV3IEFycmF5KGQrMSksIHRydWUpXG4gICAgYm91bmRhcnlbaV0gPSBjZWxsXG4gICAgbGlzdFtpXSA9IGNlbGxcbiAgfVxuICBsaXN0W2QrMV0gPSBpbml0aWFsQ2VsbFxuICBmb3IodmFyIGk9MDsgaTw9ZDsgKytpKSB7XG4gICAgdmFyIHZlcnRzID0gYm91bmRhcnlbaV0udmVydGljZXNcbiAgICB2YXIgYWRqID0gYm91bmRhcnlbaV0uYWRqYWNlbnRcbiAgICBmb3IodmFyIGo9MDsgajw9ZDsgKytqKSB7XG4gICAgICB2YXIgdiA9IHZlcnRzW2pdXG4gICAgICBpZih2IDwgMCkge1xuICAgICAgICBhZGpbal0gPSBpbml0aWFsQ2VsbFxuICAgICAgICBjb250aW51ZVxuICAgICAgfVxuICAgICAgZm9yKHZhciBrPTA7IGs8PWQ7ICsraykge1xuICAgICAgICBpZihib3VuZGFyeVtrXS52ZXJ0aWNlcy5pbmRleE9mKHYpIDwgMCkge1xuICAgICAgICAgIGFkaltqXSA9IGJvdW5kYXJ5W2tdXG4gICAgICAgIH1cbiAgICAgIH1cbiAgICB9XG4gIH1cblxuICAvL0luaXRpYWxpemUgdHJpYW5nbGVzXG4gIHZhciB0cmlhbmdsZXMgPSBuZXcgVHJpYW5ndWxhdGlvbihkLCBpbml0aWFsU2ltcGxleCwgbGlzdClcblxuICAvL0luc2VydCByZW1haW5pbmcgcG9pbnRzXG4gIHZhciB1c2VSYW5kb20gPSAhIXJhbmRvbVNlYXJjaFxuICBmb3IodmFyIGk9ZCsxOyBpPG47ICsraSkge1xuICAgIHRyaWFuZ2xlcy5pbnNlcnQocG9pbnRzW2ldLCB1c2VSYW5kb20pXG4gIH1cbiAgXG4gIC8vRXh0cmFjdCBib3VuZGFyeSBjZWxsc1xuICByZXR1cm4gdHJpYW5nbGVzLmJvdW5kYXJ5KClcbn0iLCJcInVzZSBzdHJpY3RcIlxuXG5tb2R1bGUuZXhwb3J0cyA9IGZhc3RUd29TdW1cblxuZnVuY3Rpb24gZmFzdFR3b1N1bShhLCBiLCByZXN1bHQpIHtcblx0dmFyIHggPSBhICsgYlxuXHR2YXIgYnYgPSB4IC0gYVxuXHR2YXIgYXYgPSB4IC0gYnZcblx0dmFyIGJyID0gYiAtIGJ2XG5cdHZhciBhciA9IGEgLSBhdlxuXHRpZihyZXN1bHQpIHtcblx0XHRyZXN1bHRbMF0gPSBhciArIGJyXG5cdFx0cmVzdWx0WzFdID0geFxuXHRcdHJldHVybiByZXN1bHRcblx0fVxuXHRyZXR1cm4gW2FyK2JyLCB4XVxufSIsIlwidXNlIHN0cmljdFwiXG5cbnZhciB0d29Qcm9kdWN0ID0gcmVxdWlyZShcInR3by1wcm9kdWN0XCIpXG52YXIgdHdvU3VtID0gcmVxdWlyZShcInR3by1zdW1cIilcblxubW9kdWxlLmV4cG9ydHMgPSBzY2FsZUxpbmVhckV4cGFuc2lvblxuXG5mdW5jdGlvbiBzY2FsZUxpbmVhckV4cGFuc2lvbihlLCBzY2FsZSkge1xuICB2YXIgbiA9IGUubGVuZ3RoXG4gIGlmKG4gPT09IDEpIHtcbiAgICB2YXIgdHMgPSB0d29Qcm9kdWN0KGVbMF0sIHNjYWxlKVxuICAgIGlmKHRzWzBdKSB7XG4gICAgICByZXR1cm4gdHNcbiAgICB9XG4gICAgcmV0dXJuIFsgdHNbMV0gXVxuICB9XG4gIHZhciBnID0gbmV3IEFycmF5KDIgKiBuKVxuICB2YXIgcSA9IFswLjEsIDAuMV1cbiAgdmFyIHQgPSBbMC4xLCAwLjFdXG4gIHZhciBjb3VudCA9IDBcbiAgdHdvUHJvZHVjdChlWzBdLCBzY2FsZSwgcSlcbiAgaWYocVswXSkge1xuICAgIGdbY291bnQrK10gPSBxWzBdXG4gIH1cbiAgZm9yKHZhciBpPTE7IGk8bjsgKytpKSB7XG4gICAgdHdvUHJvZHVjdChlW2ldLCBzY2FsZSwgdClcbiAgICB2YXIgcHEgPSBxWzFdXG4gICAgdHdvU3VtKHBxLCB0WzBdLCBxKVxuICAgIGlmKHFbMF0pIHtcbiAgICAgIGdbY291bnQrK10gPSBxWzBdXG4gICAgfVxuICAgIHZhciBhID0gdFsxXVxuICAgIHZhciBiID0gcVsxXVxuICAgIHZhciB4ID0gYSArIGJcbiAgICB2YXIgYnYgPSB4IC0gYVxuICAgIHZhciB5ID0gYiAtIGJ2XG4gICAgcVsxXSA9IHhcbiAgICBpZih5KSB7XG4gICAgICBnW2NvdW50KytdID0geVxuICAgIH1cbiAgfVxuICBpZihxWzFdKSB7XG4gICAgZ1tjb3VudCsrXSA9IHFbMV1cbiAgfVxuICBpZihjb3VudCA9PT0gMCkge1xuICAgIGdbY291bnQrK10gPSAwLjBcbiAgfVxuICBnLmxlbmd0aCA9IGNvdW50XG4gIHJldHVybiBnXG59IiwiXCJ1c2Ugc3RyaWN0XCJcblxubW9kdWxlLmV4cG9ydHMgPSByb2J1c3RTdWJ0cmFjdFxuXG4vL0Vhc3kgY2FzZTogQWRkIHR3byBzY2FsYXJzXG5mdW5jdGlvbiBzY2FsYXJTY2FsYXIoYSwgYikge1xuICB2YXIgeCA9IGEgKyBiXG4gIHZhciBidiA9IHggLSBhXG4gIHZhciBhdiA9IHggLSBidlxuICB2YXIgYnIgPSBiIC0gYnZcbiAgdmFyIGFyID0gYSAtIGF2XG4gIHZhciB5ID0gYXIgKyBiclxuICBpZih5KSB7XG4gICAgcmV0dXJuIFt5LCB4XVxuICB9XG4gIHJldHVybiBbeF1cbn1cblxuZnVuY3Rpb24gcm9idXN0U3VidHJhY3QoZSwgZikge1xuICB2YXIgbmUgPSBlLmxlbmd0aHwwXG4gIHZhciBuZiA9IGYubGVuZ3RofDBcbiAgaWYobmUgPT09IDEgJiYgbmYgPT09IDEpIHtcbiAgICByZXR1cm4gc2NhbGFyU2NhbGFyKGVbMF0sIC1mWzBdKVxuICB9XG4gIHZhciBuID0gbmUgKyBuZlxuICB2YXIgZyA9IG5ldyBBcnJheShuKVxuICB2YXIgY291bnQgPSAwXG4gIHZhciBlcHRyID0gMFxuICB2YXIgZnB0ciA9IDBcbiAgdmFyIGFicyA9IE1hdGguYWJzXG4gIHZhciBlaSA9IGVbZXB0cl1cbiAgdmFyIGVhID0gYWJzKGVpKVxuICB2YXIgZmkgPSAtZltmcHRyXVxuICB2YXIgZmEgPSBhYnMoZmkpXG4gIHZhciBhLCBiXG4gIGlmKGVhIDwgZmEpIHtcbiAgICBiID0gZWlcbiAgICBlcHRyICs9IDFcbiAgICBpZihlcHRyIDwgbmUpIHtcbiAgICAgIGVpID0gZVtlcHRyXVxuICAgICAgZWEgPSBhYnMoZWkpXG4gICAgfVxuICB9IGVsc2Uge1xuICAgIGIgPSBmaVxuICAgIGZwdHIgKz0gMVxuICAgIGlmKGZwdHIgPCBuZikge1xuICAgICAgZmkgPSAtZltmcHRyXVxuICAgICAgZmEgPSBhYnMoZmkpXG4gICAgfVxuICB9XG4gIGlmKChlcHRyIDwgbmUgJiYgZWEgPCBmYSkgfHwgKGZwdHIgPj0gbmYpKSB7XG4gICAgYSA9IGVpXG4gICAgZXB0ciArPSAxXG4gICAgaWYoZXB0ciA8IG5lKSB7XG4gICAgICBlaSA9IGVbZXB0cl1cbiAgICAgIGVhID0gYWJzKGVpKVxuICAgIH1cbiAgfSBlbHNlIHtcbiAgICBhID0gZmlcbiAgICBmcHRyICs9IDFcbiAgICBpZihmcHRyIDwgbmYpIHtcbiAgICAgIGZpID0gLWZbZnB0cl1cbiAgICAgIGZhID0gYWJzKGZpKVxuICAgIH1cbiAgfVxuICB2YXIgeCA9IGEgKyBiXG4gIHZhciBidiA9IHggLSBhXG4gIHZhciB5ID0gYiAtIGJ2XG4gIHZhciBxMCA9IHlcbiAgdmFyIHExID0geFxuICB2YXIgX3gsIF9idiwgX2F2LCBfYnIsIF9hclxuICB3aGlsZShlcHRyIDwgbmUgJiYgZnB0ciA8IG5mKSB7XG4gICAgaWYoZWEgPCBmYSkge1xuICAgICAgYSA9IGVpXG4gICAgICBlcHRyICs9IDFcbiAgICAgIGlmKGVwdHIgPCBuZSkge1xuICAgICAgICBlaSA9IGVbZXB0cl1cbiAgICAgICAgZWEgPSBhYnMoZWkpXG4gICAgICB9XG4gICAgfSBlbHNlIHtcbiAgICAgIGEgPSBmaVxuICAgICAgZnB0ciArPSAxXG4gICAgICBpZihmcHRyIDwgbmYpIHtcbiAgICAgICAgZmkgPSAtZltmcHRyXVxuICAgICAgICBmYSA9IGFicyhmaSlcbiAgICAgIH1cbiAgICB9XG4gICAgYiA9IHEwXG4gICAgeCA9IGEgKyBiXG4gICAgYnYgPSB4IC0gYVxuICAgIHkgPSBiIC0gYnZcbiAgICBpZih5KSB7XG4gICAgICBnW2NvdW50KytdID0geVxuICAgIH1cbiAgICBfeCA9IHExICsgeFxuICAgIF9idiA9IF94IC0gcTFcbiAgICBfYXYgPSBfeCAtIF9idlxuICAgIF9iciA9IHggLSBfYnZcbiAgICBfYXIgPSBxMSAtIF9hdlxuICAgIHEwID0gX2FyICsgX2JyXG4gICAgcTEgPSBfeFxuICB9XG4gIHdoaWxlKGVwdHIgPCBuZSkge1xuICAgIGEgPSBlaVxuICAgIGIgPSBxMFxuICAgIHggPSBhICsgYlxuICAgIGJ2ID0geCAtIGFcbiAgICB5ID0gYiAtIGJ2XG4gICAgaWYoeSkge1xuICAgICAgZ1tjb3VudCsrXSA9IHlcbiAgICB9XG4gICAgX3ggPSBxMSArIHhcbiAgICBfYnYgPSBfeCAtIHExXG4gICAgX2F2ID0gX3ggLSBfYnZcbiAgICBfYnIgPSB4IC0gX2J2XG4gICAgX2FyID0gcTEgLSBfYXZcbiAgICBxMCA9IF9hciArIF9iclxuICAgIHExID0gX3hcbiAgICBlcHRyICs9IDFcbiAgICBpZihlcHRyIDwgbmUpIHtcbiAgICAgIGVpID0gZVtlcHRyXVxuICAgIH1cbiAgfVxuICB3aGlsZShmcHRyIDwgbmYpIHtcbiAgICBhID0gZmlcbiAgICBiID0gcTBcbiAgICB4ID0gYSArIGJcbiAgICBidiA9IHggLSBhXG4gICAgeSA9IGIgLSBidlxuICAgIGlmKHkpIHtcbiAgICAgIGdbY291bnQrK10gPSB5XG4gICAgfSBcbiAgICBfeCA9IHExICsgeFxuICAgIF9idiA9IF94IC0gcTFcbiAgICBfYXYgPSBfeCAtIF9idlxuICAgIF9iciA9IHggLSBfYnZcbiAgICBfYXIgPSBxMSAtIF9hdlxuICAgIHEwID0gX2FyICsgX2JyXG4gICAgcTEgPSBfeFxuICAgIGZwdHIgKz0gMVxuICAgIGlmKGZwdHIgPCBuZikge1xuICAgICAgZmkgPSAtZltmcHRyXVxuICAgIH1cbiAgfVxuICBpZihxMCkge1xuICAgIGdbY291bnQrK10gPSBxMFxuICB9XG4gIGlmKHExKSB7XG4gICAgZ1tjb3VudCsrXSA9IHExXG4gIH1cbiAgaWYoIWNvdW50KSB7XG4gICAgZ1tjb3VudCsrXSA9IDAuMCAgXG4gIH1cbiAgZy5sZW5ndGggPSBjb3VudFxuICByZXR1cm4gZ1xufSIsIlwidXNlIHN0cmljdFwiXG5cbm1vZHVsZS5leHBvcnRzID0gbGluZWFyRXhwYW5zaW9uU3VtXG5cbi8vRWFzeSBjYXNlOiBBZGQgdHdvIHNjYWxhcnNcbmZ1bmN0aW9uIHNjYWxhclNjYWxhcihhLCBiKSB7XG4gIHZhciB4ID0gYSArIGJcbiAgdmFyIGJ2ID0geCAtIGFcbiAgdmFyIGF2ID0geCAtIGJ2XG4gIHZhciBiciA9IGIgLSBidlxuICB2YXIgYXIgPSBhIC0gYXZcbiAgdmFyIHkgPSBhciArIGJyXG4gIGlmKHkpIHtcbiAgICByZXR1cm4gW3ksIHhdXG4gIH1cbiAgcmV0dXJuIFt4XVxufVxuXG5mdW5jdGlvbiBsaW5lYXJFeHBhbnNpb25TdW0oZSwgZikge1xuICB2YXIgbmUgPSBlLmxlbmd0aHwwXG4gIHZhciBuZiA9IGYubGVuZ3RofDBcbiAgaWYobmUgPT09IDEgJiYgbmYgPT09IDEpIHtcbiAgICByZXR1cm4gc2NhbGFyU2NhbGFyKGVbMF0sIGZbMF0pXG4gIH1cbiAgdmFyIG4gPSBuZSArIG5mXG4gIHZhciBnID0gbmV3IEFycmF5KG4pXG4gIHZhciBjb3VudCA9IDBcbiAgdmFyIGVwdHIgPSAwXG4gIHZhciBmcHRyID0gMFxuICB2YXIgYWJzID0gTWF0aC5hYnNcbiAgdmFyIGVpID0gZVtlcHRyXVxuICB2YXIgZWEgPSBhYnMoZWkpXG4gIHZhciBmaSA9IGZbZnB0cl1cbiAgdmFyIGZhID0gYWJzKGZpKVxuICB2YXIgYSwgYlxuICBpZihlYSA8IGZhKSB7XG4gICAgYiA9IGVpXG4gICAgZXB0ciArPSAxXG4gICAgaWYoZXB0ciA8IG5lKSB7XG4gICAgICBlaSA9IGVbZXB0cl1cbiAgICAgIGVhID0gYWJzKGVpKVxuICAgIH1cbiAgfSBlbHNlIHtcbiAgICBiID0gZmlcbiAgICBmcHRyICs9IDFcbiAgICBpZihmcHRyIDwgbmYpIHtcbiAgICAgIGZpID0gZltmcHRyXVxuICAgICAgZmEgPSBhYnMoZmkpXG4gICAgfVxuICB9XG4gIGlmKChlcHRyIDwgbmUgJiYgZWEgPCBmYSkgfHwgKGZwdHIgPj0gbmYpKSB7XG4gICAgYSA9IGVpXG4gICAgZXB0ciArPSAxXG4gICAgaWYoZXB0ciA8IG5lKSB7XG4gICAgICBlaSA9IGVbZXB0cl1cbiAgICAgIGVhID0gYWJzKGVpKVxuICAgIH1cbiAgfSBlbHNlIHtcbiAgICBhID0gZmlcbiAgICBmcHRyICs9IDFcbiAgICBpZihmcHRyIDwgbmYpIHtcbiAgICAgIGZpID0gZltmcHRyXVxuICAgICAgZmEgPSBhYnMoZmkpXG4gICAgfVxuICB9XG4gIHZhciB4ID0gYSArIGJcbiAgdmFyIGJ2ID0geCAtIGFcbiAgdmFyIHkgPSBiIC0gYnZcbiAgdmFyIHEwID0geVxuICB2YXIgcTEgPSB4XG4gIHZhciBfeCwgX2J2LCBfYXYsIF9iciwgX2FyXG4gIHdoaWxlKGVwdHIgPCBuZSAmJiBmcHRyIDwgbmYpIHtcbiAgICBpZihlYSA8IGZhKSB7XG4gICAgICBhID0gZWlcbiAgICAgIGVwdHIgKz0gMVxuICAgICAgaWYoZXB0ciA8IG5lKSB7XG4gICAgICAgIGVpID0gZVtlcHRyXVxuICAgICAgICBlYSA9IGFicyhlaSlcbiAgICAgIH1cbiAgICB9IGVsc2Uge1xuICAgICAgYSA9IGZpXG4gICAgICBmcHRyICs9IDFcbiAgICAgIGlmKGZwdHIgPCBuZikge1xuICAgICAgICBmaSA9IGZbZnB0cl1cbiAgICAgICAgZmEgPSBhYnMoZmkpXG4gICAgICB9XG4gICAgfVxuICAgIGIgPSBxMFxuICAgIHggPSBhICsgYlxuICAgIGJ2ID0geCAtIGFcbiAgICB5ID0gYiAtIGJ2XG4gICAgaWYoeSkge1xuICAgICAgZ1tjb3VudCsrXSA9IHlcbiAgICB9XG4gICAgX3ggPSBxMSArIHhcbiAgICBfYnYgPSBfeCAtIHExXG4gICAgX2F2ID0gX3ggLSBfYnZcbiAgICBfYnIgPSB4IC0gX2J2XG4gICAgX2FyID0gcTEgLSBfYXZcbiAgICBxMCA9IF9hciArIF9iclxuICAgIHExID0gX3hcbiAgfVxuICB3aGlsZShlcHRyIDwgbmUpIHtcbiAgICBhID0gZWlcbiAgICBiID0gcTBcbiAgICB4ID0gYSArIGJcbiAgICBidiA9IHggLSBhXG4gICAgeSA9IGIgLSBidlxuICAgIGlmKHkpIHtcbiAgICAgIGdbY291bnQrK10gPSB5XG4gICAgfVxuICAgIF94ID0gcTEgKyB4XG4gICAgX2J2ID0gX3ggLSBxMVxuICAgIF9hdiA9IF94IC0gX2J2XG4gICAgX2JyID0geCAtIF9idlxuICAgIF9hciA9IHExIC0gX2F2XG4gICAgcTAgPSBfYXIgKyBfYnJcbiAgICBxMSA9IF94XG4gICAgZXB0ciArPSAxXG4gICAgaWYoZXB0ciA8IG5lKSB7XG4gICAgICBlaSA9IGVbZXB0cl1cbiAgICB9XG4gIH1cbiAgd2hpbGUoZnB0ciA8IG5mKSB7XG4gICAgYSA9IGZpXG4gICAgYiA9IHEwXG4gICAgeCA9IGEgKyBiXG4gICAgYnYgPSB4IC0gYVxuICAgIHkgPSBiIC0gYnZcbiAgICBpZih5KSB7XG4gICAgICBnW2NvdW50KytdID0geVxuICAgIH0gXG4gICAgX3ggPSBxMSArIHhcbiAgICBfYnYgPSBfeCAtIHExXG4gICAgX2F2ID0gX3ggLSBfYnZcbiAgICBfYnIgPSB4IC0gX2J2XG4gICAgX2FyID0gcTEgLSBfYXZcbiAgICBxMCA9IF9hciArIF9iclxuICAgIHExID0gX3hcbiAgICBmcHRyICs9IDFcbiAgICBpZihmcHRyIDwgbmYpIHtcbiAgICAgIGZpID0gZltmcHRyXVxuICAgIH1cbiAgfVxuICBpZihxMCkge1xuICAgIGdbY291bnQrK10gPSBxMFxuICB9XG4gIGlmKHExKSB7XG4gICAgZ1tjb3VudCsrXSA9IHExXG4gIH1cbiAgaWYoIWNvdW50KSB7XG4gICAgZ1tjb3VudCsrXSA9IDAuMCAgXG4gIH1cbiAgZy5sZW5ndGggPSBjb3VudFxuICByZXR1cm4gZ1xufSIsIlwidXNlIHN0cmljdFwiXG5cbm1vZHVsZS5leHBvcnRzID0gdHdvUHJvZHVjdFxuXG52YXIgU1BMSVRURVIgPSArKE1hdGgucG93KDIsIDI3KSArIDEuMClcblxuZnVuY3Rpb24gdHdvUHJvZHVjdChhLCBiLCByZXN1bHQpIHtcbiAgdmFyIHggPSBhICogYlxuXG4gIHZhciBjID0gU1BMSVRURVIgKiBhXG4gIHZhciBhYmlnID0gYyAtIGFcbiAgdmFyIGFoaSA9IGMgLSBhYmlnXG4gIHZhciBhbG8gPSBhIC0gYWhpXG5cbiAgdmFyIGQgPSBTUExJVFRFUiAqIGJcbiAgdmFyIGJiaWcgPSBkIC0gYlxuICB2YXIgYmhpID0gZCAtIGJiaWdcbiAgdmFyIGJsbyA9IGIgLSBiaGlcblxuICB2YXIgZXJyMSA9IHggLSAoYWhpICogYmhpKVxuICB2YXIgZXJyMiA9IGVycjEgLSAoYWxvICogYmhpKVxuICB2YXIgZXJyMyA9IGVycjIgLSAoYWhpICogYmxvKVxuXG4gIHZhciB5ID0gYWxvICogYmxvIC0gZXJyM1xuXG4gIGlmKHJlc3VsdCkge1xuICAgIHJlc3VsdFswXSA9IHlcbiAgICByZXN1bHRbMV0gPSB4XG4gICAgcmV0dXJuIHJlc3VsdFxuICB9XG5cbiAgcmV0dXJuIFsgeSwgeCBdXG59IiwiXCJ1c2Ugc3RyaWN0XCJcblxudmFyIHR3b1Byb2R1Y3QgPSByZXF1aXJlKFwidHdvLXByb2R1Y3RcIilcbnZhciByb2J1c3RTdW0gPSByZXF1aXJlKFwicm9idXN0LXN1bVwiKVxudmFyIHJvYnVzdFNjYWxlID0gcmVxdWlyZShcInJvYnVzdC1zY2FsZVwiKVxudmFyIHJvYnVzdFN1YnRyYWN0ID0gcmVxdWlyZShcInJvYnVzdC1zdWJ0cmFjdFwiKVxuXG52YXIgTlVNX0VYUEFORCA9IDVcblxudmFyIEVQU0lMT04gICAgID0gMS4xMTAyMjMwMjQ2MjUxNTY1ZS0xNlxudmFyIEVSUkJPVU5EMyAgID0gKDMuMCArIDE2LjAgKiBFUFNJTE9OKSAqIEVQU0lMT05cbnZhciBFUlJCT1VORDQgICA9ICg3LjAgKyA1Ni4wICogRVBTSUxPTikgKiBFUFNJTE9OXG5cbmZ1bmN0aW9uIGNvZmFjdG9yKG0sIGMpIHtcbiAgdmFyIHJlc3VsdCA9IG5ldyBBcnJheShtLmxlbmd0aC0xKVxuICBmb3IodmFyIGk9MTsgaTxtLmxlbmd0aDsgKytpKSB7XG4gICAgdmFyIHIgPSByZXN1bHRbaS0xXSA9IG5ldyBBcnJheShtLmxlbmd0aC0xKVxuICAgIGZvcih2YXIgaj0wLGs9MDsgajxtLmxlbmd0aDsgKytqKSB7XG4gICAgICBpZihqID09PSBjKSB7XG4gICAgICAgIGNvbnRpbnVlXG4gICAgICB9XG4gICAgICByW2srK10gPSBtW2ldW2pdXG4gICAgfVxuICB9XG4gIHJldHVybiByZXN1bHRcbn1cblxuZnVuY3Rpb24gbWF0cml4KG4pIHtcbiAgdmFyIHJlc3VsdCA9IG5ldyBBcnJheShuKVxuICBmb3IodmFyIGk9MDsgaTxuOyArK2kpIHtcbiAgICByZXN1bHRbaV0gPSBuZXcgQXJyYXkobilcbiAgICBmb3IodmFyIGo9MDsgajxuOyArK2opIHtcbiAgICAgIHJlc3VsdFtpXVtqXSA9IFtcIm1cIiwgaiwgXCJbXCIsIChuLWktMSksIFwiXVwiXS5qb2luKFwiXCIpXG4gICAgfVxuICB9XG4gIHJldHVybiByZXN1bHRcbn1cblxuZnVuY3Rpb24gc2lnbihuKSB7XG4gIGlmKG4gJiAxKSB7XG4gICAgcmV0dXJuIFwiLVwiXG4gIH1cbiAgcmV0dXJuIFwiXCJcbn1cblxuZnVuY3Rpb24gZ2VuZXJhdGVTdW0oZXhwcikge1xuICBpZihleHByLmxlbmd0aCA9PT0gMSkge1xuICAgIHJldHVybiBleHByWzBdXG4gIH0gZWxzZSBpZihleHByLmxlbmd0aCA9PT0gMikge1xuICAgIHJldHVybiBbXCJzdW0oXCIsIGV4cHJbMF0sIFwiLFwiLCBleHByWzFdLCBcIilcIl0uam9pbihcIlwiKVxuICB9IGVsc2Uge1xuICAgIHZhciBtID0gZXhwci5sZW5ndGg+PjFcbiAgICByZXR1cm4gW1wic3VtKFwiLCBnZW5lcmF0ZVN1bShleHByLnNsaWNlKDAsIG0pKSwgXCIsXCIsIGdlbmVyYXRlU3VtKGV4cHIuc2xpY2UobSkpLCBcIilcIl0uam9pbihcIlwiKVxuICB9XG59XG5cbmZ1bmN0aW9uIGRldGVybWluYW50KG0pIHtcbiAgaWYobS5sZW5ndGggPT09IDIpIHtcbiAgICByZXR1cm4gW1tcInN1bShwcm9kKFwiLCBtWzBdWzBdLCBcIixcIiwgbVsxXVsxXSwgXCIpLHByb2QoLVwiLCBtWzBdWzFdLCBcIixcIiwgbVsxXVswXSwgXCIpKVwiXS5qb2luKFwiXCIpXVxuICB9IGVsc2Uge1xuICAgIHZhciBleHByID0gW11cbiAgICBmb3IodmFyIGk9MDsgaTxtLmxlbmd0aDsgKytpKSB7XG4gICAgICBleHByLnB1c2goW1wic2NhbGUoXCIsIGdlbmVyYXRlU3VtKGRldGVybWluYW50KGNvZmFjdG9yKG0sIGkpKSksIFwiLFwiLCBzaWduKGkpLCBtWzBdW2ldLCBcIilcIl0uam9pbihcIlwiKSlcbiAgICB9XG4gICAgcmV0dXJuIGV4cHJcbiAgfVxufVxuXG5mdW5jdGlvbiBvcmllbnRhdGlvbihuKSB7XG4gIHZhciBwb3MgPSBbXVxuICB2YXIgbmVnID0gW11cbiAgdmFyIG0gPSBtYXRyaXgobilcbiAgdmFyIGFyZ3MgPSBbXVxuICBmb3IodmFyIGk9MDsgaTxuOyArK2kpIHtcbiAgICBpZigoaSYxKT09PTApIHtcbiAgICAgIHBvcy5wdXNoLmFwcGx5KHBvcywgZGV0ZXJtaW5hbnQoY29mYWN0b3IobSwgaSkpKVxuICAgIH0gZWxzZSB7XG4gICAgICBuZWcucHVzaC5hcHBseShuZWcsIGRldGVybWluYW50KGNvZmFjdG9yKG0sIGkpKSlcbiAgICB9XG4gICAgYXJncy5wdXNoKFwibVwiICsgaSlcbiAgfVxuICB2YXIgcG9zRXhwciA9IGdlbmVyYXRlU3VtKHBvcylcbiAgdmFyIG5lZ0V4cHIgPSBnZW5lcmF0ZVN1bShuZWcpXG4gIHZhciBmdW5jTmFtZSA9IFwib3JpZW50YXRpb25cIiArIG4gKyBcIkV4YWN0XCJcbiAgdmFyIGNvZGUgPSBbXCJmdW5jdGlvbiBcIiwgZnVuY05hbWUsIFwiKFwiLCBhcmdzLmpvaW4oKSwgXCIpe3ZhciBwPVwiLCBwb3NFeHByLCBcIixuPVwiLCBuZWdFeHByLCBcIixkPXN1YihwLG4pO1xcXG5yZXR1cm4gZFtkLmxlbmd0aC0xXTt9O3JldHVybiBcIiwgZnVuY05hbWVdLmpvaW4oXCJcIilcbiAgdmFyIHByb2MgPSBuZXcgRnVuY3Rpb24oXCJzdW1cIiwgXCJwcm9kXCIsIFwic2NhbGVcIiwgXCJzdWJcIiwgY29kZSlcbiAgcmV0dXJuIHByb2Mocm9idXN0U3VtLCB0d29Qcm9kdWN0LCByb2J1c3RTY2FsZSwgcm9idXN0U3VidHJhY3QpXG59XG5cbnZhciBvcmllbnRhdGlvbjNFeGFjdCA9IG9yaWVudGF0aW9uKDMpXG52YXIgb3JpZW50YXRpb240RXhhY3QgPSBvcmllbnRhdGlvbig0KVxuXG52YXIgQ0FDSEVEID0gW1xuICBmdW5jdGlvbiBvcmllbnRhdGlvbjAoKSB7IHJldHVybiAwIH0sXG4gIGZ1bmN0aW9uIG9yaWVudGF0aW9uMSgpIHsgcmV0dXJuIDAgfSxcbiAgZnVuY3Rpb24gb3JpZW50YXRpb24yKGEsIGIpIHsgXG4gICAgcmV0dXJuIGJbMF0gLSBhWzBdXG4gIH0sXG4gIGZ1bmN0aW9uIG9yaWVudGF0aW9uMyhhLCBiLCBjKSB7XG4gICAgdmFyIGwgPSAoYVsxXSAtIGNbMV0pICogKGJbMF0gLSBjWzBdKVxuICAgIHZhciByID0gKGFbMF0gLSBjWzBdKSAqIChiWzFdIC0gY1sxXSlcbiAgICB2YXIgZGV0ID0gbCAtIHJcbiAgICB2YXIgc1xuICAgIGlmKGwgPiAwKSB7XG4gICAgICBpZihyIDw9IDApIHtcbiAgICAgICAgcmV0dXJuIGRldFxuICAgICAgfSBlbHNlIHtcbiAgICAgICAgcyA9IGwgKyByXG4gICAgICB9XG4gICAgfSBlbHNlIGlmKGwgPCAwKSB7XG4gICAgICBpZihyID49IDApIHtcbiAgICAgICAgcmV0dXJuIGRldFxuICAgICAgfSBlbHNlIHtcbiAgICAgICAgcyA9IC0obCArIHIpXG4gICAgICB9XG4gICAgfSBlbHNlIHtcbiAgICAgIHJldHVybiBkZXRcbiAgICB9XG4gICAgdmFyIHRvbCA9IEVSUkJPVU5EMyAqIHNcbiAgICBpZihkZXQgPj0gdG9sIHx8IGRldCA8PSAtdG9sKSB7XG4gICAgICByZXR1cm4gZGV0XG4gICAgfVxuICAgIHJldHVybiBvcmllbnRhdGlvbjNFeGFjdChhLCBiLCBjKVxuICB9LFxuICBmdW5jdGlvbiBvcmllbnRhdGlvbjQoYSxiLGMsZCkge1xuICAgIHZhciBhZHggPSBhWzBdIC0gZFswXVxuICAgIHZhciBiZHggPSBiWzBdIC0gZFswXVxuICAgIHZhciBjZHggPSBjWzBdIC0gZFswXVxuICAgIHZhciBhZHkgPSBhWzFdIC0gZFsxXVxuICAgIHZhciBiZHkgPSBiWzFdIC0gZFsxXVxuICAgIHZhciBjZHkgPSBjWzFdIC0gZFsxXVxuICAgIHZhciBhZHogPSBhWzJdIC0gZFsyXVxuICAgIHZhciBiZHogPSBiWzJdIC0gZFsyXVxuICAgIHZhciBjZHogPSBjWzJdIC0gZFsyXVxuICAgIHZhciBiZHhjZHkgPSBiZHggKiBjZHlcbiAgICB2YXIgY2R4YmR5ID0gY2R4ICogYmR5XG4gICAgdmFyIGNkeGFkeSA9IGNkeCAqIGFkeVxuICAgIHZhciBhZHhjZHkgPSBhZHggKiBjZHlcbiAgICB2YXIgYWR4YmR5ID0gYWR4ICogYmR5XG4gICAgdmFyIGJkeGFkeSA9IGJkeCAqIGFkeVxuICAgIHZhciBkZXQgPSBhZHogKiAoYmR4Y2R5IC0gY2R4YmR5KSBcbiAgICAgICAgICAgICsgYmR6ICogKGNkeGFkeSAtIGFkeGNkeSlcbiAgICAgICAgICAgICsgY2R6ICogKGFkeGJkeSAtIGJkeGFkeSlcbiAgICB2YXIgcGVybWFuZW50ID0gKE1hdGguYWJzKGJkeGNkeSkgKyBNYXRoLmFicyhjZHhiZHkpKSAqIE1hdGguYWJzKGFkeilcbiAgICAgICAgICAgICAgICAgICsgKE1hdGguYWJzKGNkeGFkeSkgKyBNYXRoLmFicyhhZHhjZHkpKSAqIE1hdGguYWJzKGJkeilcbiAgICAgICAgICAgICAgICAgICsgKE1hdGguYWJzKGFkeGJkeSkgKyBNYXRoLmFicyhiZHhhZHkpKSAqIE1hdGguYWJzKGNkeilcbiAgICB2YXIgdG9sID0gRVJSQk9VTkQ0ICogcGVybWFuZW50XG4gICAgaWYgKChkZXQgPiB0b2wpIHx8ICgtZGV0ID4gdG9sKSkge1xuICAgICAgcmV0dXJuIGRldFxuICAgIH1cbiAgICByZXR1cm4gb3JpZW50YXRpb240RXhhY3QoYSxiLGMsZClcbiAgfVxuXVxuXG5mdW5jdGlvbiBzbG93T3JpZW50KGFyZ3MpIHtcbiAgdmFyIHByb2MgPSBDQUNIRURbYXJncy5sZW5ndGhdXG4gIGlmKCFwcm9jKSB7XG4gICAgcHJvYyA9IENBQ0hFRFthcmdzLmxlbmd0aF0gPSBvcmllbnRhdGlvbihhcmdzLmxlbmd0aClcbiAgfVxuICByZXR1cm4gcHJvYy5hcHBseSh1bmRlZmluZWQsIGFyZ3MpXG59XG5cbmZ1bmN0aW9uIGdlbmVyYXRlT3JpZW50YXRpb25Qcm9jKCkge1xuICB3aGlsZShDQUNIRUQubGVuZ3RoIDw9IE5VTV9FWFBBTkQpIHtcbiAgICBDQUNIRUQucHVzaChvcmllbnRhdGlvbihDQUNIRUQubGVuZ3RoKSlcbiAgfVxuICB2YXIgYXJncyA9IFtdXG4gIHZhciBwcm9jQXJncyA9IFtcInNsb3dcIl1cbiAgZm9yKHZhciBpPTA7IGk8PU5VTV9FWFBBTkQ7ICsraSkge1xuICAgIGFyZ3MucHVzaChcImFcIiArIGkpXG4gICAgcHJvY0FyZ3MucHVzaChcIm9cIiArIGkpXG4gIH1cbiAgdmFyIGNvZGUgPSBbXG4gICAgXCJmdW5jdGlvbiBnZXRPcmllbnRhdGlvbihcIiwgYXJncy5qb2luKCksIFwiKXtzd2l0Y2goYXJndW1lbnRzLmxlbmd0aCl7Y2FzZSAwOmNhc2UgMTpyZXR1cm4gMDtcIlxuICBdXG4gIGZvcih2YXIgaT0yOyBpPD1OVU1fRVhQQU5EOyArK2kpIHtcbiAgICBjb2RlLnB1c2goXCJjYXNlIFwiLCBpLCBcIjpyZXR1cm4gb1wiLCBpLCBcIihcIiwgYXJncy5zbGljZSgwLCBpKS5qb2luKCksIFwiKTtcIilcbiAgfVxuICBjb2RlLnB1c2goXCJ9dmFyIHM9bmV3IEFycmF5KGFyZ3VtZW50cy5sZW5ndGgpO2Zvcih2YXIgaT0wO2k8YXJndW1lbnRzLmxlbmd0aDsrK2kpe3NbaV09YXJndW1lbnRzW2ldfTtyZXR1cm4gc2xvdyhzKTt9cmV0dXJuIGdldE9yaWVudGF0aW9uXCIpXG4gIHByb2NBcmdzLnB1c2goY29kZS5qb2luKFwiXCIpKVxuXG4gIHZhciBwcm9jID0gRnVuY3Rpb24uYXBwbHkodW5kZWZpbmVkLCBwcm9jQXJncylcbiAgbW9kdWxlLmV4cG9ydHMgPSBwcm9jLmFwcGx5KHVuZGVmaW5lZCwgW3Nsb3dPcmllbnRdLmNvbmNhdChDQUNIRUQpKVxuICBmb3IodmFyIGk9MDsgaTw9TlVNX0VYUEFORDsgKytpKSB7XG4gICAgbW9kdWxlLmV4cG9ydHNbaV0gPSBDQUNIRURbaV1cbiAgfVxufVxuXG5nZW5lcmF0ZU9yaWVudGF0aW9uUHJvYygpIiwiLyoqXG4gKiBCaXQgdHdpZGRsaW5nIGhhY2tzIGZvciBKYXZhU2NyaXB0LlxuICpcbiAqIEF1dGhvcjogTWlrb2xhIEx5c2Vua29cbiAqXG4gKiBQb3J0ZWQgZnJvbSBTdGFuZm9yZCBiaXQgdHdpZGRsaW5nIGhhY2sgbGlicmFyeTpcbiAqICAgIGh0dHA6Ly9ncmFwaGljcy5zdGFuZm9yZC5lZHUvfnNlYW5kZXIvYml0aGFja3MuaHRtbFxuICovXG5cblwidXNlIHN0cmljdFwiOyBcInVzZSByZXN0cmljdFwiO1xuXG4vL051bWJlciBvZiBiaXRzIGluIGFuIGludGVnZXJcbnZhciBJTlRfQklUUyA9IDMyO1xuXG4vL0NvbnN0YW50c1xuZXhwb3J0cy5JTlRfQklUUyAgPSBJTlRfQklUUztcbmV4cG9ydHMuSU5UX01BWCAgID0gIDB4N2ZmZmZmZmY7XG5leHBvcnRzLklOVF9NSU4gICA9IC0xPDwoSU5UX0JJVFMtMSk7XG5cbi8vUmV0dXJucyAtMSwgMCwgKzEgZGVwZW5kaW5nIG9uIHNpZ24gb2YgeFxuZXhwb3J0cy5zaWduID0gZnVuY3Rpb24odikge1xuICByZXR1cm4gKHYgPiAwKSAtICh2IDwgMCk7XG59XG5cbi8vQ29tcHV0ZXMgYWJzb2x1dGUgdmFsdWUgb2YgaW50ZWdlclxuZXhwb3J0cy5hYnMgPSBmdW5jdGlvbih2KSB7XG4gIHZhciBtYXNrID0gdiA+PiAoSU5UX0JJVFMtMSk7XG4gIHJldHVybiAodiBeIG1hc2spIC0gbWFzaztcbn1cblxuLy9Db21wdXRlcyBtaW5pbXVtIG9mIGludGVnZXJzIHggYW5kIHlcbmV4cG9ydHMubWluID0gZnVuY3Rpb24oeCwgeSkge1xuICByZXR1cm4geSBeICgoeCBeIHkpICYgLSh4IDwgeSkpO1xufVxuXG4vL0NvbXB1dGVzIG1heGltdW0gb2YgaW50ZWdlcnMgeCBhbmQgeVxuZXhwb3J0cy5tYXggPSBmdW5jdGlvbih4LCB5KSB7XG4gIHJldHVybiB4IF4gKCh4IF4geSkgJiAtKHggPCB5KSk7XG59XG5cbi8vQ2hlY2tzIGlmIGEgbnVtYmVyIGlzIGEgcG93ZXIgb2YgdHdvXG5leHBvcnRzLmlzUG93MiA9IGZ1bmN0aW9uKHYpIHtcbiAgcmV0dXJuICEodiAmICh2LTEpKSAmJiAoISF2KTtcbn1cblxuLy9Db21wdXRlcyBsb2cgYmFzZSAyIG9mIHZcbmV4cG9ydHMubG9nMiA9IGZ1bmN0aW9uKHYpIHtcbiAgdmFyIHIsIHNoaWZ0O1xuICByID0gICAgICh2ID4gMHhGRkZGKSA8PCA0OyB2ID4+Pj0gcjtcbiAgc2hpZnQgPSAodiA+IDB4RkYgICkgPDwgMzsgdiA+Pj49IHNoaWZ0OyByIHw9IHNoaWZ0O1xuICBzaGlmdCA9ICh2ID4gMHhGICAgKSA8PCAyOyB2ID4+Pj0gc2hpZnQ7IHIgfD0gc2hpZnQ7XG4gIHNoaWZ0ID0gKHYgPiAweDMgICApIDw8IDE7IHYgPj4+PSBzaGlmdDsgciB8PSBzaGlmdDtcbiAgcmV0dXJuIHIgfCAodiA+PiAxKTtcbn1cblxuLy9Db21wdXRlcyBsb2cgYmFzZSAxMCBvZiB2XG5leHBvcnRzLmxvZzEwID0gZnVuY3Rpb24odikge1xuICByZXR1cm4gICh2ID49IDEwMDAwMDAwMDApID8gOSA6ICh2ID49IDEwMDAwMDAwMCkgPyA4IDogKHYgPj0gMTAwMDAwMDApID8gNyA6XG4gICAgICAgICAgKHYgPj0gMTAwMDAwMCkgPyA2IDogKHYgPj0gMTAwMDAwKSA/IDUgOiAodiA+PSAxMDAwMCkgPyA0IDpcbiAgICAgICAgICAodiA+PSAxMDAwKSA/IDMgOiAodiA+PSAxMDApID8gMiA6ICh2ID49IDEwKSA/IDEgOiAwO1xufVxuXG4vL0NvdW50cyBudW1iZXIgb2YgYml0c1xuZXhwb3J0cy5wb3BDb3VudCA9IGZ1bmN0aW9uKHYpIHtcbiAgdiA9IHYgLSAoKHYgPj4+IDEpICYgMHg1NTU1NTU1NSk7XG4gIHYgPSAodiAmIDB4MzMzMzMzMzMpICsgKCh2ID4+PiAyKSAmIDB4MzMzMzMzMzMpO1xuICByZXR1cm4gKCh2ICsgKHYgPj4+IDQpICYgMHhGMEYwRjBGKSAqIDB4MTAxMDEwMSkgPj4+IDI0O1xufVxuXG4vL0NvdW50cyBudW1iZXIgb2YgdHJhaWxpbmcgemVyb3NcbmZ1bmN0aW9uIGNvdW50VHJhaWxpbmdaZXJvcyh2KSB7XG4gIHZhciBjID0gMzI7XG4gIHYgJj0gLXY7XG4gIGlmICh2KSBjLS07XG4gIGlmICh2ICYgMHgwMDAwRkZGRikgYyAtPSAxNjtcbiAgaWYgKHYgJiAweDAwRkYwMEZGKSBjIC09IDg7XG4gIGlmICh2ICYgMHgwRjBGMEYwRikgYyAtPSA0O1xuICBpZiAodiAmIDB4MzMzMzMzMzMpIGMgLT0gMjtcbiAgaWYgKHYgJiAweDU1NTU1NTU1KSBjIC09IDE7XG4gIHJldHVybiBjO1xufVxuZXhwb3J0cy5jb3VudFRyYWlsaW5nWmVyb3MgPSBjb3VudFRyYWlsaW5nWmVyb3M7XG5cbi8vUm91bmRzIHRvIG5leHQgcG93ZXIgb2YgMlxuZXhwb3J0cy5uZXh0UG93MiA9IGZ1bmN0aW9uKHYpIHtcbiAgdiArPSB2ID09PSAwO1xuICAtLXY7XG4gIHYgfD0gdiA+Pj4gMTtcbiAgdiB8PSB2ID4+PiAyO1xuICB2IHw9IHYgPj4+IDQ7XG4gIHYgfD0gdiA+Pj4gODtcbiAgdiB8PSB2ID4+PiAxNjtcbiAgcmV0dXJuIHYgKyAxO1xufVxuXG4vL1JvdW5kcyBkb3duIHRvIHByZXZpb3VzIHBvd2VyIG9mIDJcbmV4cG9ydHMucHJldlBvdzIgPSBmdW5jdGlvbih2KSB7XG4gIHYgfD0gdiA+Pj4gMTtcbiAgdiB8PSB2ID4+PiAyO1xuICB2IHw9IHYgPj4+IDQ7XG4gIHYgfD0gdiA+Pj4gODtcbiAgdiB8PSB2ID4+PiAxNjtcbiAgcmV0dXJuIHYgLSAodj4+PjEpO1xufVxuXG4vL0NvbXB1dGVzIHBhcml0eSBvZiB3b3JkXG5leHBvcnRzLnBhcml0eSA9IGZ1bmN0aW9uKHYpIHtcbiAgdiBePSB2ID4+PiAxNjtcbiAgdiBePSB2ID4+PiA4O1xuICB2IF49IHYgPj4+IDQ7XG4gIHYgJj0gMHhmO1xuICByZXR1cm4gKDB4Njk5NiA+Pj4gdikgJiAxO1xufVxuXG52YXIgUkVWRVJTRV9UQUJMRSA9IG5ldyBBcnJheSgyNTYpO1xuXG4oZnVuY3Rpb24odGFiKSB7XG4gIGZvcih2YXIgaT0wOyBpPDI1NjsgKytpKSB7XG4gICAgdmFyIHYgPSBpLCByID0gaSwgcyA9IDc7XG4gICAgZm9yICh2ID4+Pj0gMTsgdjsgdiA+Pj49IDEpIHtcbiAgICAgIHIgPDw9IDE7XG4gICAgICByIHw9IHYgJiAxO1xuICAgICAgLS1zO1xuICAgIH1cbiAgICB0YWJbaV0gPSAociA8PCBzKSAmIDB4ZmY7XG4gIH1cbn0pKFJFVkVSU0VfVEFCTEUpO1xuXG4vL1JldmVyc2UgYml0cyBpbiBhIDMyIGJpdCB3b3JkXG5leHBvcnRzLnJldmVyc2UgPSBmdW5jdGlvbih2KSB7XG4gIHJldHVybiAgKFJFVkVSU0VfVEFCTEVbIHYgICAgICAgICAmIDB4ZmZdIDw8IDI0KSB8XG4gICAgICAgICAgKFJFVkVSU0VfVEFCTEVbKHYgPj4+IDgpICAmIDB4ZmZdIDw8IDE2KSB8XG4gICAgICAgICAgKFJFVkVSU0VfVEFCTEVbKHYgPj4+IDE2KSAmIDB4ZmZdIDw8IDgpICB8XG4gICAgICAgICAgIFJFVkVSU0VfVEFCTEVbKHYgPj4+IDI0KSAmIDB4ZmZdO1xufVxuXG4vL0ludGVybGVhdmUgYml0cyBvZiAyIGNvb3JkaW5hdGVzIHdpdGggMTYgYml0cy4gIFVzZWZ1bCBmb3IgZmFzdCBxdWFkdHJlZSBjb2Rlc1xuZXhwb3J0cy5pbnRlcmxlYXZlMiA9IGZ1bmN0aW9uKHgsIHkpIHtcbiAgeCAmPSAweEZGRkY7XG4gIHggPSAoeCB8ICh4IDw8IDgpKSAmIDB4MDBGRjAwRkY7XG4gIHggPSAoeCB8ICh4IDw8IDQpKSAmIDB4MEYwRjBGMEY7XG4gIHggPSAoeCB8ICh4IDw8IDIpKSAmIDB4MzMzMzMzMzM7XG4gIHggPSAoeCB8ICh4IDw8IDEpKSAmIDB4NTU1NTU1NTU7XG5cbiAgeSAmPSAweEZGRkY7XG4gIHkgPSAoeSB8ICh5IDw8IDgpKSAmIDB4MDBGRjAwRkY7XG4gIHkgPSAoeSB8ICh5IDw8IDQpKSAmIDB4MEYwRjBGMEY7XG4gIHkgPSAoeSB8ICh5IDw8IDIpKSAmIDB4MzMzMzMzMzM7XG4gIHkgPSAoeSB8ICh5IDw8IDEpKSAmIDB4NTU1NTU1NTU7XG5cbiAgcmV0dXJuIHggfCAoeSA8PCAxKTtcbn1cblxuLy9FeHRyYWN0cyB0aGUgbnRoIGludGVybGVhdmVkIGNvbXBvbmVudFxuZXhwb3J0cy5kZWludGVybGVhdmUyID0gZnVuY3Rpb24odiwgbikge1xuICB2ID0gKHYgPj4+IG4pICYgMHg1NTU1NTU1NTtcbiAgdiA9ICh2IHwgKHYgPj4+IDEpKSAgJiAweDMzMzMzMzMzO1xuICB2ID0gKHYgfCAodiA+Pj4gMikpICAmIDB4MEYwRjBGMEY7XG4gIHYgPSAodiB8ICh2ID4+PiA0KSkgICYgMHgwMEZGMDBGRjtcbiAgdiA9ICh2IHwgKHYgPj4+IDE2KSkgJiAweDAwMEZGRkY7XG4gIHJldHVybiAodiA8PCAxNikgPj4gMTY7XG59XG5cblxuLy9JbnRlcmxlYXZlIGJpdHMgb2YgMyBjb29yZGluYXRlcywgZWFjaCB3aXRoIDEwIGJpdHMuICBVc2VmdWwgZm9yIGZhc3Qgb2N0cmVlIGNvZGVzXG5leHBvcnRzLmludGVybGVhdmUzID0gZnVuY3Rpb24oeCwgeSwgeikge1xuICB4ICY9IDB4M0ZGO1xuICB4ICA9ICh4IHwgKHg8PDE2KSkgJiA0Mjc4MTkwMzM1O1xuICB4ICA9ICh4IHwgKHg8PDgpKSAgJiAyNTE3MTk2OTU7XG4gIHggID0gKHggfCAoeDw8NCkpICAmIDMyNzIzNTYwMzU7XG4gIHggID0gKHggfCAoeDw8MikpICAmIDEyMjcxMzM1MTM7XG5cbiAgeSAmPSAweDNGRjtcbiAgeSAgPSAoeSB8ICh5PDwxNikpICYgNDI3ODE5MDMzNTtcbiAgeSAgPSAoeSB8ICh5PDw4KSkgICYgMjUxNzE5Njk1O1xuICB5ICA9ICh5IHwgKHk8PDQpKSAgJiAzMjcyMzU2MDM1O1xuICB5ICA9ICh5IHwgKHk8PDIpKSAgJiAxMjI3MTMzNTEzO1xuICB4IHw9ICh5IDw8IDEpO1xuICBcbiAgeiAmPSAweDNGRjtcbiAgeiAgPSAoeiB8ICh6PDwxNikpICYgNDI3ODE5MDMzNTtcbiAgeiAgPSAoeiB8ICh6PDw4KSkgICYgMjUxNzE5Njk1O1xuICB6ICA9ICh6IHwgKHo8PDQpKSAgJiAzMjcyMzU2MDM1O1xuICB6ICA9ICh6IHwgKHo8PDIpKSAgJiAxMjI3MTMzNTEzO1xuICBcbiAgcmV0dXJuIHggfCAoeiA8PCAyKTtcbn1cblxuLy9FeHRyYWN0cyBudGggaW50ZXJsZWF2ZWQgY29tcG9uZW50IG9mIGEgMy10dXBsZVxuZXhwb3J0cy5kZWludGVybGVhdmUzID0gZnVuY3Rpb24odiwgbikge1xuICB2ID0gKHYgPj4+IG4pICAgICAgICYgMTIyNzEzMzUxMztcbiAgdiA9ICh2IHwgKHY+Pj4yKSkgICAmIDMyNzIzNTYwMzU7XG4gIHYgPSAodiB8ICh2Pj4+NCkpICAgJiAyNTE3MTk2OTU7XG4gIHYgPSAodiB8ICh2Pj4+OCkpICAgJiA0Mjc4MTkwMzM1O1xuICB2ID0gKHYgfCAodj4+PjE2KSkgICYgMHgzRkY7XG4gIHJldHVybiAodjw8MjIpPj4yMjtcbn1cblxuLy9Db21wdXRlcyBuZXh0IGNvbWJpbmF0aW9uIGluIGNvbGV4aWNvZ3JhcGhpYyBvcmRlciAodGhpcyBpcyBtaXN0YWtlbmx5IGNhbGxlZCBuZXh0UGVybXV0YXRpb24gb24gdGhlIGJpdCB0d2lkZGxpbmcgaGFja3MgcGFnZSlcbmV4cG9ydHMubmV4dENvbWJpbmF0aW9uID0gZnVuY3Rpb24odikge1xuICB2YXIgdCA9IHYgfCAodiAtIDEpO1xuICByZXR1cm4gKHQgKyAxKSB8ICgoKH50ICYgLX50KSAtIDEpID4+PiAoY291bnRUcmFpbGluZ1plcm9zKHYpICsgMSkpO1xufVxuXG4iLCJcInVzZSBzdHJpY3RcIjsgXCJ1c2UgcmVzdHJpY3RcIjtcblxubW9kdWxlLmV4cG9ydHMgPSBVbmlvbkZpbmQ7XG5cbmZ1bmN0aW9uIFVuaW9uRmluZChjb3VudCkge1xuICB0aGlzLnJvb3RzID0gbmV3IEFycmF5KGNvdW50KTtcbiAgdGhpcy5yYW5rcyA9IG5ldyBBcnJheShjb3VudCk7XG4gIFxuICBmb3IodmFyIGk9MDsgaTxjb3VudDsgKytpKSB7XG4gICAgdGhpcy5yb290c1tpXSA9IGk7XG4gICAgdGhpcy5yYW5rc1tpXSA9IDA7XG4gIH1cbn1cblxudmFyIHByb3RvID0gVW5pb25GaW5kLnByb3RvdHlwZVxuXG5PYmplY3QuZGVmaW5lUHJvcGVydHkocHJvdG8sIFwibGVuZ3RoXCIsIHtcbiAgXCJnZXRcIjogZnVuY3Rpb24oKSB7XG4gICAgcmV0dXJuIHRoaXMucm9vdHMubGVuZ3RoXG4gIH1cbn0pXG5cbnByb3RvLm1ha2VTZXQgPSBmdW5jdGlvbigpIHtcbiAgdmFyIG4gPSB0aGlzLnJvb3RzLmxlbmd0aDtcbiAgdGhpcy5yb290cy5wdXNoKG4pO1xuICB0aGlzLnJhbmtzLnB1c2goMCk7XG4gIHJldHVybiBuO1xufVxuXG5wcm90by5maW5kID0gZnVuY3Rpb24oeCkge1xuICB2YXIgcm9vdHMgPSB0aGlzLnJvb3RzO1xuICB3aGlsZShyb290c1t4XSAhPT0geCkge1xuICAgIHZhciB5ID0gcm9vdHNbeF07XG4gICAgcm9vdHNbeF0gPSByb290c1t5XTtcbiAgICB4ID0geTtcbiAgfVxuICByZXR1cm4geDtcbn1cblxucHJvdG8ubGluayA9IGZ1bmN0aW9uKHgsIHkpIHtcbiAgdmFyIHhyID0gdGhpcy5maW5kKHgpXG4gICAgLCB5ciA9IHRoaXMuZmluZCh5KTtcbiAgaWYoeHIgPT09IHlyKSB7XG4gICAgcmV0dXJuO1xuICB9XG4gIHZhciByYW5rcyA9IHRoaXMucmFua3NcbiAgICAsIHJvb3RzID0gdGhpcy5yb290c1xuICAgICwgeGQgICAgPSByYW5rc1t4cl1cbiAgICAsIHlkICAgID0gcmFua3NbeXJdO1xuICBpZih4ZCA8IHlkKSB7XG4gICAgcm9vdHNbeHJdID0geXI7XG4gIH0gZWxzZSBpZih5ZCA8IHhkKSB7XG4gICAgcm9vdHNbeXJdID0geHI7XG4gIH0gZWxzZSB7XG4gICAgcm9vdHNbeXJdID0geHI7XG4gICAgKytyYW5rc1t4cl07XG4gIH1cbn0iLCJcInVzZSBzdHJpY3RcIjsgXCJ1c2UgcmVzdHJpY3RcIjtcblxudmFyIGJpdHMgICAgICA9IHJlcXVpcmUoXCJiaXQtdHdpZGRsZVwiKVxuICAsIFVuaW9uRmluZCA9IHJlcXVpcmUoXCJ1bmlvbi1maW5kXCIpXG5cbi8vUmV0dXJucyB0aGUgZGltZW5zaW9uIG9mIGEgY2VsbCBjb21wbGV4XG5mdW5jdGlvbiBkaW1lbnNpb24oY2VsbHMpIHtcbiAgdmFyIGQgPSAwXG4gICAgLCBtYXggPSBNYXRoLm1heFxuICBmb3IodmFyIGk9MCwgaWw9Y2VsbHMubGVuZ3RoOyBpPGlsOyArK2kpIHtcbiAgICBkID0gbWF4KGQsIGNlbGxzW2ldLmxlbmd0aClcbiAgfVxuICByZXR1cm4gZC0xXG59XG5leHBvcnRzLmRpbWVuc2lvbiA9IGRpbWVuc2lvblxuXG4vL0NvdW50cyB0aGUgbnVtYmVyIG9mIHZlcnRpY2VzIGluIGZhY2VzXG5mdW5jdGlvbiBjb3VudFZlcnRpY2VzKGNlbGxzKSB7XG4gIHZhciB2YyA9IC0xXG4gICAgLCBtYXggPSBNYXRoLm1heFxuICBmb3IodmFyIGk9MCwgaWw9Y2VsbHMubGVuZ3RoOyBpPGlsOyArK2kpIHtcbiAgICB2YXIgYyA9IGNlbGxzW2ldXG4gICAgZm9yKHZhciBqPTAsIGpsPWMubGVuZ3RoOyBqPGpsOyArK2opIHtcbiAgICAgIHZjID0gbWF4KHZjLCBjW2pdKVxuICAgIH1cbiAgfVxuICByZXR1cm4gdmMrMVxufVxuZXhwb3J0cy5jb3VudFZlcnRpY2VzID0gY291bnRWZXJ0aWNlc1xuXG4vL1JldHVybnMgYSBkZWVwIGNvcHkgb2YgY2VsbHNcbmZ1bmN0aW9uIGNsb25lQ2VsbHMoY2VsbHMpIHtcbiAgdmFyIG5jZWxscyA9IG5ldyBBcnJheShjZWxscy5sZW5ndGgpXG4gIGZvcih2YXIgaT0wLCBpbD1jZWxscy5sZW5ndGg7IGk8aWw7ICsraSkge1xuICAgIG5jZWxsc1tpXSA9IGNlbGxzW2ldLnNsaWNlKDApXG4gIH1cbiAgcmV0dXJuIG5jZWxsc1xufVxuZXhwb3J0cy5jbG9uZUNlbGxzID0gY2xvbmVDZWxsc1xuXG4vL1JhbmtzIGEgcGFpciBvZiBjZWxscyB1cCB0byBwZXJtdXRhdGlvblxuZnVuY3Rpb24gY29tcGFyZUNlbGxzKGEsIGIpIHtcbiAgdmFyIG4gPSBhLmxlbmd0aFxuICAgICwgdCA9IGEubGVuZ3RoIC0gYi5sZW5ndGhcbiAgICAsIG1pbiA9IE1hdGgubWluXG4gIGlmKHQpIHtcbiAgICByZXR1cm4gdFxuICB9XG4gIHN3aXRjaChuKSB7XG4gICAgY2FzZSAwOlxuICAgICAgcmV0dXJuIDA7XG4gICAgY2FzZSAxOlxuICAgICAgcmV0dXJuIGFbMF0gLSBiWzBdO1xuICAgIGNhc2UgMjpcbiAgICAgIHZhciBkID0gYVswXSthWzFdLWJbMF0tYlsxXVxuICAgICAgaWYoZCkge1xuICAgICAgICByZXR1cm4gZFxuICAgICAgfVxuICAgICAgcmV0dXJuIG1pbihhWzBdLGFbMV0pIC0gbWluKGJbMF0sYlsxXSlcbiAgICBjYXNlIDM6XG4gICAgICB2YXIgbDEgPSBhWzBdK2FbMV1cbiAgICAgICAgLCBtMSA9IGJbMF0rYlsxXVxuICAgICAgZCA9IGwxK2FbMl0gLSAobTErYlsyXSlcbiAgICAgIGlmKGQpIHtcbiAgICAgICAgcmV0dXJuIGRcbiAgICAgIH1cbiAgICAgIHZhciBsMCA9IG1pbihhWzBdLCBhWzFdKVxuICAgICAgICAsIG0wID0gbWluKGJbMF0sIGJbMV0pXG4gICAgICAgICwgZCAgPSBtaW4obDAsIGFbMl0pIC0gbWluKG0wLCBiWzJdKVxuICAgICAgaWYoZCkge1xuICAgICAgICByZXR1cm4gZFxuICAgICAgfVxuICAgICAgcmV0dXJuIG1pbihsMCthWzJdLCBsMSkgLSBtaW4obTArYlsyXSwgbTEpXG4gICAgXG4gICAgLy9UT0RPOiBNYXliZSBvcHRpbWl6ZSBuPTQgYXMgd2VsbD9cbiAgICBcbiAgICBkZWZhdWx0OlxuICAgICAgdmFyIGFzID0gYS5zbGljZSgwKVxuICAgICAgYXMuc29ydCgpXG4gICAgICB2YXIgYnMgPSBiLnNsaWNlKDApXG4gICAgICBicy5zb3J0KClcbiAgICAgIGZvcih2YXIgaT0wOyBpPG47ICsraSkge1xuICAgICAgICB0ID0gYXNbaV0gLSBic1tpXVxuICAgICAgICBpZih0KSB7XG4gICAgICAgICAgcmV0dXJuIHRcbiAgICAgICAgfVxuICAgICAgfVxuICAgICAgcmV0dXJuIDBcbiAgfVxufVxuZXhwb3J0cy5jb21wYXJlQ2VsbHMgPSBjb21wYXJlQ2VsbHNcblxuZnVuY3Rpb24gY29tcGFyZVppcHBlZChhLCBiKSB7XG4gIHJldHVybiBjb21wYXJlQ2VsbHMoYVswXSwgYlswXSlcbn1cblxuLy9QdXRzIGEgY2VsbCBjb21wbGV4IGludG8gbm9ybWFsIG9yZGVyIGZvciB0aGUgcHVycG9zZXMgb2YgZmluZENlbGwgcXVlcmllc1xuZnVuY3Rpb24gbm9ybWFsaXplKGNlbGxzLCBhdHRyKSB7XG4gIGlmKGF0dHIpIHtcbiAgICB2YXIgbGVuID0gY2VsbHMubGVuZ3RoXG4gICAgdmFyIHppcHBlZCA9IG5ldyBBcnJheShsZW4pXG4gICAgZm9yKHZhciBpPTA7IGk8bGVuOyArK2kpIHtcbiAgICAgIHppcHBlZFtpXSA9IFtjZWxsc1tpXSwgYXR0cltpXV1cbiAgICB9XG4gICAgemlwcGVkLnNvcnQoY29tcGFyZVppcHBlZClcbiAgICBmb3IodmFyIGk9MDsgaTxsZW47ICsraSkge1xuICAgICAgY2VsbHNbaV0gPSB6aXBwZWRbaV1bMF1cbiAgICAgIGF0dHJbaV0gPSB6aXBwZWRbaV1bMV1cbiAgICB9XG4gICAgcmV0dXJuIGNlbGxzXG4gIH0gZWxzZSB7XG4gICAgY2VsbHMuc29ydChjb21wYXJlQ2VsbHMpXG4gICAgcmV0dXJuIGNlbGxzXG4gIH1cbn1cbmV4cG9ydHMubm9ybWFsaXplID0gbm9ybWFsaXplXG5cbi8vUmVtb3ZlcyBhbGwgZHVwbGljYXRlIGNlbGxzIGluIHRoZSBjb21wbGV4XG5mdW5jdGlvbiB1bmlxdWUoY2VsbHMpIHtcbiAgaWYoY2VsbHMubGVuZ3RoID09PSAwKSB7XG4gICAgcmV0dXJuIFtdXG4gIH1cbiAgdmFyIHB0ciA9IDFcbiAgICAsIGxlbiA9IGNlbGxzLmxlbmd0aFxuICBmb3IodmFyIGk9MTsgaTxsZW47ICsraSkge1xuICAgIHZhciBhID0gY2VsbHNbaV1cbiAgICBpZihjb21wYXJlQ2VsbHMoYSwgY2VsbHNbaS0xXSkpIHtcbiAgICAgIGlmKGkgPT09IHB0cikge1xuICAgICAgICBwdHIrK1xuICAgICAgICBjb250aW51ZVxuICAgICAgfVxuICAgICAgY2VsbHNbcHRyKytdID0gYVxuICAgIH1cbiAgfVxuICBjZWxscy5sZW5ndGggPSBwdHJcbiAgcmV0dXJuIGNlbGxzXG59XG5leHBvcnRzLnVuaXF1ZSA9IHVuaXF1ZTtcblxuLy9GaW5kcyBhIGNlbGwgaW4gYSBub3JtYWxpemVkIGNlbGwgY29tcGxleFxuZnVuY3Rpb24gZmluZENlbGwoY2VsbHMsIGMpIHtcbiAgdmFyIGxvID0gMFxuICAgICwgaGkgPSBjZWxscy5sZW5ndGgtMVxuICAgICwgciAgPSAtMVxuICB3aGlsZSAobG8gPD0gaGkpIHtcbiAgICB2YXIgbWlkID0gKGxvICsgaGkpID4+IDFcbiAgICAgICwgcyAgID0gY29tcGFyZUNlbGxzKGNlbGxzW21pZF0sIGMpXG4gICAgaWYocyA8PSAwKSB7XG4gICAgICBpZihzID09PSAwKSB7XG4gICAgICAgIHIgPSBtaWRcbiAgICAgIH1cbiAgICAgIGxvID0gbWlkICsgMVxuICAgIH0gZWxzZSBpZihzID4gMCkge1xuICAgICAgaGkgPSBtaWQgLSAxXG4gICAgfVxuICB9XG4gIHJldHVybiByXG59XG5leHBvcnRzLmZpbmRDZWxsID0gZmluZENlbGw7XG5cbi8vQnVpbGRzIGFuIGluZGV4IGZvciBhbiBuLWNlbGwuICBUaGlzIGlzIG1vcmUgZ2VuZXJhbCB0aGFuIGR1YWwsIGJ1dCBsZXNzIGVmZmljaWVudFxuZnVuY3Rpb24gaW5jaWRlbmNlKGZyb21fY2VsbHMsIHRvX2NlbGxzKSB7XG4gIHZhciBpbmRleCA9IG5ldyBBcnJheShmcm9tX2NlbGxzLmxlbmd0aClcbiAgZm9yKHZhciBpPTAsIGlsPWluZGV4Lmxlbmd0aDsgaTxpbDsgKytpKSB7XG4gICAgaW5kZXhbaV0gPSBbXVxuICB9XG4gIHZhciBiID0gW11cbiAgZm9yKHZhciBpPTAsIG49dG9fY2VsbHMubGVuZ3RoOyBpPG47ICsraSkge1xuICAgIHZhciBjID0gdG9fY2VsbHNbaV1cbiAgICB2YXIgY2wgPSBjLmxlbmd0aFxuICAgIGZvcih2YXIgaz0xLCBrbj0oMTw8Y2wpOyBrPGtuOyArK2spIHtcbiAgICAgIGIubGVuZ3RoID0gYml0cy5wb3BDb3VudChrKVxuICAgICAgdmFyIGwgPSAwXG4gICAgICBmb3IodmFyIGo9MDsgajxjbDsgKytqKSB7XG4gICAgICAgIGlmKGsgJiAoMTw8aikpIHtcbiAgICAgICAgICBiW2wrK10gPSBjW2pdXG4gICAgICAgIH1cbiAgICAgIH1cbiAgICAgIHZhciBpZHg9ZmluZENlbGwoZnJvbV9jZWxscywgYilcbiAgICAgIGlmKGlkeCA8IDApIHtcbiAgICAgICAgY29udGludWVcbiAgICAgIH1cbiAgICAgIHdoaWxlKHRydWUpIHtcbiAgICAgICAgaW5kZXhbaWR4KytdLnB1c2goaSlcbiAgICAgICAgaWYoaWR4ID49IGZyb21fY2VsbHMubGVuZ3RoIHx8IGNvbXBhcmVDZWxscyhmcm9tX2NlbGxzW2lkeF0sIGIpICE9PSAwKSB7XG4gICAgICAgICAgYnJlYWtcbiAgICAgICAgfVxuICAgICAgfVxuICAgIH1cbiAgfVxuICByZXR1cm4gaW5kZXhcbn1cbmV4cG9ydHMuaW5jaWRlbmNlID0gaW5jaWRlbmNlXG5cbi8vQ29tcHV0ZXMgdGhlIGR1YWwgb2YgdGhlIG1lc2guICBUaGlzIGlzIGJhc2ljYWxseSBhbiBvcHRpbWl6ZWQgdmVyc2lvbiBvZiBidWlsZEluZGV4IGZvciB0aGUgc2l0dWF0aW9uIHdoZXJlIGZyb21fY2VsbHMgaXMganVzdCB0aGUgbGlzdCBvZiB2ZXJ0aWNlc1xuZnVuY3Rpb24gZHVhbChjZWxscywgdmVydGV4X2NvdW50KSB7XG4gIGlmKCF2ZXJ0ZXhfY291bnQpIHtcbiAgICByZXR1cm4gaW5jaWRlbmNlKHVuaXF1ZShza2VsZXRvbihjZWxscywgMCkpLCBjZWxscywgMClcbiAgfVxuICB2YXIgcmVzID0gbmV3IEFycmF5KHZlcnRleF9jb3VudClcbiAgZm9yKHZhciBpPTA7IGk8dmVydGV4X2NvdW50OyArK2kpIHtcbiAgICByZXNbaV0gPSBbXVxuICB9XG4gIGZvcih2YXIgaT0wLCBsZW49Y2VsbHMubGVuZ3RoOyBpPGxlbjsgKytpKSB7XG4gICAgdmFyIGMgPSBjZWxsc1tpXVxuICAgIGZvcih2YXIgaj0wLCBjbD1jLmxlbmd0aDsgajxjbDsgKytqKSB7XG4gICAgICByZXNbY1tqXV0ucHVzaChpKVxuICAgIH1cbiAgfVxuICByZXR1cm4gcmVzXG59XG5leHBvcnRzLmR1YWwgPSBkdWFsXG5cbi8vRW51bWVyYXRlcyBhbGwgY2VsbHMgaW4gdGhlIGNvbXBsZXhcbmZ1bmN0aW9uIGV4cGxvZGUoY2VsbHMpIHtcbiAgdmFyIHJlc3VsdCA9IFtdXG4gIGZvcih2YXIgaT0wLCBpbD1jZWxscy5sZW5ndGg7IGk8aWw7ICsraSkge1xuICAgIHZhciBjID0gY2VsbHNbaV1cbiAgICAgICwgY2wgPSBjLmxlbmd0aHwwXG4gICAgZm9yKHZhciBqPTEsIGpsPSgxPDxjbCk7IGo8amw7ICsraikge1xuICAgICAgdmFyIGIgPSBbXVxuICAgICAgZm9yKHZhciBrPTA7IGs8Y2w7ICsraykge1xuICAgICAgICBpZigoaiA+Pj4gaykgJiAxKSB7XG4gICAgICAgICAgYi5wdXNoKGNba10pXG4gICAgICAgIH1cbiAgICAgIH1cbiAgICAgIHJlc3VsdC5wdXNoKGIpXG4gICAgfVxuICB9XG4gIHJldHVybiBub3JtYWxpemUocmVzdWx0KVxufVxuZXhwb3J0cy5leHBsb2RlID0gZXhwbG9kZVxuXG4vL0VudW1lcmF0ZXMgYWxsIG9mIHRoZSBuLWNlbGxzIG9mIGEgY2VsbCBjb21wbGV4XG5mdW5jdGlvbiBza2VsZXRvbihjZWxscywgbikge1xuICBpZihuIDwgMCkge1xuICAgIHJldHVybiBbXVxuICB9XG4gIHZhciByZXN1bHQgPSBbXVxuICAgICwgazAgICAgID0gKDE8PChuKzEpKS0xXG4gIGZvcih2YXIgaT0wOyBpPGNlbGxzLmxlbmd0aDsgKytpKSB7XG4gICAgdmFyIGMgPSBjZWxsc1tpXVxuICAgIGZvcih2YXIgaz1rMDsgazwoMTw8Yy5sZW5ndGgpOyBrPWJpdHMubmV4dENvbWJpbmF0aW9uKGspKSB7XG4gICAgICB2YXIgYiA9IG5ldyBBcnJheShuKzEpXG4gICAgICAgICwgbCA9IDBcbiAgICAgIGZvcih2YXIgaj0wOyBqPGMubGVuZ3RoOyArK2opIHtcbiAgICAgICAgaWYoayAmICgxPDxqKSkge1xuICAgICAgICAgIGJbbCsrXSA9IGNbal1cbiAgICAgICAgfVxuICAgICAgfVxuICAgICAgcmVzdWx0LnB1c2goYilcbiAgICB9XG4gIH1cbiAgcmV0dXJuIG5vcm1hbGl6ZShyZXN1bHQpXG59XG5leHBvcnRzLnNrZWxldG9uID0gc2tlbGV0b247XG5cbi8vQ29tcHV0ZXMgdGhlIGJvdW5kYXJ5IG9mIGFsbCBjZWxscywgZG9lcyBub3QgcmVtb3ZlIGR1cGxpY2F0ZXNcbmZ1bmN0aW9uIGJvdW5kYXJ5KGNlbGxzKSB7XG4gIHZhciByZXMgPSBbXVxuICBmb3IodmFyIGk9MCxpbD1jZWxscy5sZW5ndGg7IGk8aWw7ICsraSkge1xuICAgIHZhciBjID0gY2VsbHNbaV1cbiAgICBmb3IodmFyIGo9MCxjbD1jLmxlbmd0aDsgajxjbDsgKytqKSB7XG4gICAgICB2YXIgYiA9IG5ldyBBcnJheShjLmxlbmd0aC0xKVxuICAgICAgZm9yKHZhciBrPTAsIGw9MDsgazxjbDsgKytrKSB7XG4gICAgICAgIGlmKGsgIT09IGopIHtcbiAgICAgICAgICBiW2wrK10gPSBjW2tdXG4gICAgICAgIH1cbiAgICAgIH1cbiAgICAgIHJlcy5wdXNoKGIpXG4gICAgfVxuICB9XG4gIHJldHVybiBub3JtYWxpemUocmVzKVxufVxuZXhwb3J0cy5ib3VuZGFyeSA9IGJvdW5kYXJ5O1xuXG4vL0NvbXB1dGVzIGNvbm5lY3RlZCBjb21wb25lbnRzIGZvciBhIGRlbnNlIGNlbGwgY29tcGxleFxuZnVuY3Rpb24gY29ubmVjdGVkQ29tcG9uZW50c19kZW5zZShjZWxscywgdmVydGV4X2NvdW50KSB7XG4gIHZhciBsYWJlbHMgPSBuZXcgVW5pb25GaW5kKHZlcnRleF9jb3VudClcbiAgZm9yKHZhciBpPTA7IGk8Y2VsbHMubGVuZ3RoOyArK2kpIHtcbiAgICB2YXIgYyA9IGNlbGxzW2ldXG4gICAgZm9yKHZhciBqPTA7IGo8Yy5sZW5ndGg7ICsraikge1xuICAgICAgZm9yKHZhciBrPWorMTsgazxjLmxlbmd0aDsgKytrKSB7XG4gICAgICAgIGxhYmVscy5saW5rKGNbal0sIGNba10pXG4gICAgICB9XG4gICAgfVxuICB9XG4gIHZhciBjb21wb25lbnRzID0gW11cbiAgICAsIGNvbXBvbmVudF9sYWJlbHMgPSBsYWJlbHMucmFua3NcbiAgZm9yKHZhciBpPTA7IGk8Y29tcG9uZW50X2xhYmVscy5sZW5ndGg7ICsraSkge1xuICAgIGNvbXBvbmVudF9sYWJlbHNbaV0gPSAtMVxuICB9XG4gIGZvcih2YXIgaT0wOyBpPGNlbGxzLmxlbmd0aDsgKytpKSB7XG4gICAgdmFyIGwgPSBsYWJlbHMuZmluZChjZWxsc1tpXVswXSlcbiAgICBpZihjb21wb25lbnRfbGFiZWxzW2xdIDwgMCkge1xuICAgICAgY29tcG9uZW50X2xhYmVsc1tsXSA9IGNvbXBvbmVudHMubGVuZ3RoXG4gICAgICBjb21wb25lbnRzLnB1c2goW2NlbGxzW2ldLnNsaWNlKDApXSlcbiAgICB9IGVsc2Uge1xuICAgICAgY29tcG9uZW50c1tjb21wb25lbnRfbGFiZWxzW2xdXS5wdXNoKGNlbGxzW2ldLnNsaWNlKDApKVxuICAgIH1cbiAgfVxuICByZXR1cm4gY29tcG9uZW50c1xufVxuXG4vL0NvbXB1dGVzIGNvbm5lY3RlZCBjb21wb25lbnRzIGZvciBhIHNwYXJzZSBncmFwaFxuZnVuY3Rpb24gY29ubmVjdGVkQ29tcG9uZW50c19zcGFyc2UoY2VsbHMpIHtcbiAgdmFyIHZlcnRpY2VzICA9IHVuaXF1ZShub3JtYWxpemUoc2tlbGV0b24oY2VsbHMsIDApKSlcbiAgICAsIGxhYmVscyAgICA9IG5ldyBVbmlvbkZpbmQodmVydGljZXMubGVuZ3RoKVxuICBmb3IodmFyIGk9MDsgaTxjZWxscy5sZW5ndGg7ICsraSkge1xuICAgIHZhciBjID0gY2VsbHNbaV1cbiAgICBmb3IodmFyIGo9MDsgajxjLmxlbmd0aDsgKytqKSB7XG4gICAgICB2YXIgdmogPSBmaW5kQ2VsbCh2ZXJ0aWNlcywgW2Nbal1dKVxuICAgICAgZm9yKHZhciBrPWorMTsgazxjLmxlbmd0aDsgKytrKSB7XG4gICAgICAgIGxhYmVscy5saW5rKHZqLCBmaW5kQ2VsbCh2ZXJ0aWNlcywgW2Nba11dKSlcbiAgICAgIH1cbiAgICB9XG4gIH1cbiAgdmFyIGNvbXBvbmVudHMgICAgICAgID0gW11cbiAgICAsIGNvbXBvbmVudF9sYWJlbHMgID0gbGFiZWxzLnJhbmtzXG4gIGZvcih2YXIgaT0wOyBpPGNvbXBvbmVudF9sYWJlbHMubGVuZ3RoOyArK2kpIHtcbiAgICBjb21wb25lbnRfbGFiZWxzW2ldID0gLTFcbiAgfVxuICBmb3IodmFyIGk9MDsgaTxjZWxscy5sZW5ndGg7ICsraSkge1xuICAgIHZhciBsID0gbGFiZWxzLmZpbmQoZmluZENlbGwodmVydGljZXMsIFtjZWxsc1tpXVswXV0pKTtcbiAgICBpZihjb21wb25lbnRfbGFiZWxzW2xdIDwgMCkge1xuICAgICAgY29tcG9uZW50X2xhYmVsc1tsXSA9IGNvbXBvbmVudHMubGVuZ3RoXG4gICAgICBjb21wb25lbnRzLnB1c2goW2NlbGxzW2ldLnNsaWNlKDApXSlcbiAgICB9IGVsc2Uge1xuICAgICAgY29tcG9uZW50c1tjb21wb25lbnRfbGFiZWxzW2xdXS5wdXNoKGNlbGxzW2ldLnNsaWNlKDApKVxuICAgIH1cbiAgfVxuICByZXR1cm4gY29tcG9uZW50c1xufVxuXG4vL0NvbXB1dGVzIGNvbm5lY3RlZCBjb21wb25lbnRzIGZvciBhIGNlbGwgY29tcGxleFxuZnVuY3Rpb24gY29ubmVjdGVkQ29tcG9uZW50cyhjZWxscywgdmVydGV4X2NvdW50KSB7XG4gIGlmKHZlcnRleF9jb3VudCkge1xuICAgIHJldHVybiBjb25uZWN0ZWRDb21wb25lbnRzX2RlbnNlKGNlbGxzLCB2ZXJ0ZXhfY291bnQpXG4gIH1cbiAgcmV0dXJuIGNvbm5lY3RlZENvbXBvbmVudHNfc3BhcnNlKGNlbGxzKVxufVxuZXhwb3J0cy5jb25uZWN0ZWRDb21wb25lbnRzID0gY29ubmVjdGVkQ29tcG9uZW50c1xuIiwiXCJ1c2Ugc3RyaWN0XCJcblxuZnVuY3Rpb24gdW5pcXVlX3ByZWQobGlzdCwgY29tcGFyZSkge1xuICB2YXIgcHRyID0gMVxuICAgICwgbGVuID0gbGlzdC5sZW5ndGhcbiAgICAsIGE9bGlzdFswXSwgYj1saXN0WzBdXG4gIGZvcih2YXIgaT0xOyBpPGxlbjsgKytpKSB7XG4gICAgYiA9IGFcbiAgICBhID0gbGlzdFtpXVxuICAgIGlmKGNvbXBhcmUoYSwgYikpIHtcbiAgICAgIGlmKGkgPT09IHB0cikge1xuICAgICAgICBwdHIrK1xuICAgICAgICBjb250aW51ZVxuICAgICAgfVxuICAgICAgbGlzdFtwdHIrK10gPSBhXG4gICAgfVxuICB9XG4gIGxpc3QubGVuZ3RoID0gcHRyXG4gIHJldHVybiBsaXN0XG59XG5cbmZ1bmN0aW9uIHVuaXF1ZV9lcShsaXN0KSB7XG4gIHZhciBwdHIgPSAxXG4gICAgLCBsZW4gPSBsaXN0Lmxlbmd0aFxuICAgICwgYT1saXN0WzBdLCBiID0gbGlzdFswXVxuICBmb3IodmFyIGk9MTsgaTxsZW47ICsraSwgYj1hKSB7XG4gICAgYiA9IGFcbiAgICBhID0gbGlzdFtpXVxuICAgIGlmKGEgIT09IGIpIHtcbiAgICAgIGlmKGkgPT09IHB0cikge1xuICAgICAgICBwdHIrK1xuICAgICAgICBjb250aW51ZVxuICAgICAgfVxuICAgICAgbGlzdFtwdHIrK10gPSBhXG4gICAgfVxuICB9XG4gIGxpc3QubGVuZ3RoID0gcHRyXG4gIHJldHVybiBsaXN0XG59XG5cbmZ1bmN0aW9uIHVuaXF1ZShsaXN0LCBjb21wYXJlLCBzb3J0ZWQpIHtcbiAgaWYobGlzdC5sZW5ndGggPT09IDApIHtcbiAgICByZXR1cm4gbGlzdFxuICB9XG4gIGlmKGNvbXBhcmUpIHtcbiAgICBpZighc29ydGVkKSB7XG4gICAgICBsaXN0LnNvcnQoY29tcGFyZSlcbiAgICB9XG4gICAgcmV0dXJuIHVuaXF1ZV9wcmVkKGxpc3QsIGNvbXBhcmUpXG4gIH1cbiAgaWYoIXNvcnRlZCkge1xuICAgIGxpc3Quc29ydCgpXG4gIH1cbiAgcmV0dXJuIHVuaXF1ZV9lcShsaXN0KVxufVxuXG5tb2R1bGUuZXhwb3J0cyA9IHVuaXF1ZVxuIiwiXCJ1c2Ugc3RyaWN0XCJcblxudmFyIGNoID0gcmVxdWlyZShcImluY3JlbWVudGFsLWNvbnZleC1odWxsXCIpXG52YXIgdW5pcSA9IHJlcXVpcmUoXCJ1bmlxXCIpXG5cbm1vZHVsZS5leHBvcnRzID0gdHJpYW5ndWxhdGVcblxuZnVuY3Rpb24gTGlmdGVkUG9pbnQocCwgaSkge1xuICB0aGlzLnBvaW50ID0gcFxuICB0aGlzLmluZGV4ID0gaVxufVxuXG5mdW5jdGlvbiBjb21wYXJlTGlmdGVkKGEsIGIpIHtcbiAgdmFyIGFwID0gYS5wb2ludFxuICB2YXIgYnAgPSBiLnBvaW50XG4gIHZhciBkID0gYXAubGVuZ3RoXG4gIGZvcih2YXIgaT0wOyBpPGQ7ICsraSkge1xuICAgIHZhciBzID0gYnBbaV0gLSBhcFtpXVxuICAgIGlmKHMpIHtcbiAgICAgIHJldHVybiBzXG4gICAgfVxuICB9XG4gIHJldHVybiAwXG59XG5cbmZ1bmN0aW9uIHRyaWFuZ3VsYXRlMUQobiwgcG9pbnRzLCBpbmNsdWRlUG9pbnRBdEluZmluaXR5KSB7XG4gIGlmKG4gPT09IDEpIHtcbiAgICBpZihpbmNsdWRlUG9pbnRBdEluZmluaXR5KSB7XG4gICAgICByZXR1cm4gWyBbLTEsIDBdIF1cbiAgICB9IGVsc2Uge1xuICAgICAgcmV0dXJuIFtdXG4gICAgfVxuICB9XG4gIHZhciBsaWZ0ZWQgPSBwb2ludHMubWFwKGZ1bmN0aW9uKHAsIGkpIHtcbiAgICByZXR1cm4gWyBwWzBdLCBpIF1cbiAgfSlcbiAgbGlmdGVkLnNvcnQoZnVuY3Rpb24oYSxiKSB7XG4gICAgcmV0dXJuIGFbMF0gLSBiWzBdXG4gIH0pXG4gIHZhciBjZWxscyA9IG5ldyBBcnJheShuIC0gMSlcbiAgZm9yKHZhciBpPTE7IGk8bjsgKytpKSB7XG4gICAgdmFyIGEgPSBsaWZ0ZWRbaS0xXVxuICAgIHZhciBiID0gbGlmdGVkW2ldXG4gICAgY2VsbHNbaS0xXSA9IFsgYVsxXSwgYlsxXSBdXG4gIH1cbiAgaWYoaW5jbHVkZVBvaW50QXRJbmZpbml0eSkge1xuICAgIGNlbGxzLnB1c2goXG4gICAgICBbIC0xLCBjZWxsc1swXVsxXSwgXSxcbiAgICAgIFsgY2VsbHNbbi0xXVsxXSwgLTEgXSlcbiAgfVxuICByZXR1cm4gY2VsbHNcbn1cblxuZnVuY3Rpb24gdHJpYW5ndWxhdGUocG9pbnRzLCBpbmNsdWRlUG9pbnRBdEluZmluaXR5KSB7XG4gIHZhciBuID0gcG9pbnRzLmxlbmd0aFxuICBpZihuID09PSAwKSB7XG4gICAgcmV0dXJuIFtdXG4gIH1cbiAgXG4gIHZhciBkID0gcG9pbnRzWzBdLmxlbmd0aFxuICBpZihkIDwgMSkge1xuICAgIHJldHVybiBbXVxuICB9XG5cbiAgLy9TcGVjaWFsIGNhc2U6ICBGb3IgMUQgd2UgY2FuIGp1c3Qgc29ydCB0aGUgcG9pbnRzXG4gIGlmKGQgPT09IDEpIHtcbiAgICByZXR1cm4gdHJpYW5ndWxhdGUxRChuLCBwb2ludHMsIGluY2x1ZGVQb2ludEF0SW5maW5pdHkpXG4gIH1cbiAgXG4gIC8vTGlmdCBwb2ludHMsIHNvcnRcbiAgdmFyIGxpZnRlZCA9IG5ldyBBcnJheShuKVxuICB2YXIgdXBwZXIgPSAxLjBcbiAgZm9yKHZhciBpPTA7IGk8bjsgKytpKSB7XG4gICAgdmFyIHAgPSBwb2ludHNbaV1cbiAgICB2YXIgeCA9IG5ldyBBcnJheShkKzEpXG4gICAgdmFyIGwgPSAwLjBcbiAgICBmb3IodmFyIGo9MDsgajxkOyArK2opIHtcbiAgICAgIHZhciB2ID0gcFtqXVxuICAgICAgeFtqXSA9IHZcbiAgICAgIGwgKz0gdiAqIHZcbiAgICB9XG4gICAgeFtkXSA9IGxcbiAgICBsaWZ0ZWRbaV0gPSBuZXcgTGlmdGVkUG9pbnQoeCwgaSlcbiAgICB1cHBlciA9IE1hdGgubWF4KGwsIHVwcGVyKVxuICB9XG4gIHVuaXEobGlmdGVkLCBjb21wYXJlTGlmdGVkKVxuICBcbiAgLy9Eb3VibGUgcG9pbnRzXG4gIG4gPSBsaWZ0ZWQubGVuZ3RoXG5cbiAgLy9DcmVhdGUgbmV3IGxpc3Qgb2YgcG9pbnRzXG4gIHZhciBkcG9pbnRzID0gbmV3IEFycmF5KG4gKyBkICsgMSlcbiAgdmFyIGRpbmRleCA9IG5ldyBBcnJheShuICsgZCArIDEpXG5cbiAgLy9BZGQgc3RlaW5lciBwb2ludHMgYXQgdG9wXG4gIHZhciB1ID0gKGQrMSkgKiAoZCsxKSAqIHVwcGVyXG4gIHZhciB5ID0gbmV3IEFycmF5KGQrMSlcbiAgZm9yKHZhciBpPTA7IGk8PWQ7ICsraSkge1xuICAgIHlbaV0gPSAwLjBcbiAgfVxuICB5W2RdID0gdVxuXG4gIGRwb2ludHNbMF0gPSB5LnNsaWNlKClcbiAgZGluZGV4WzBdID0gLTFcblxuICBmb3IodmFyIGk9MDsgaTw9ZDsgKytpKSB7XG4gICAgdmFyIHggPSB5LnNsaWNlKClcbiAgICB4W2ldID0gMVxuICAgIGRwb2ludHNbaSsxXSA9IHhcbiAgICBkaW5kZXhbaSsxXSA9IC0xXG4gIH1cblxuICAvL0NvcHkgcmVzdCBvZiB0aGUgcG9pbnRzIG92ZXJcbiAgZm9yKHZhciBpPTA7IGk8bjsgKytpKSB7XG4gICAgdmFyIGggPSBsaWZ0ZWRbaV1cbiAgICBkcG9pbnRzW2kgKyBkICsgMV0gPSBoLnBvaW50XG4gICAgZGluZGV4W2kgKyBkICsgMV0gPSAgaC5pbmRleFxuICB9XG5cbiAgLy9Db25zdHJ1Y3QgY29udmV4IGh1bGxcbiAgdmFyIGh1bGwgPSBjaChkcG9pbnRzLCBmYWxzZSlcbiAgaWYoaW5jbHVkZVBvaW50QXRJbmZpbml0eSkge1xuICAgIGh1bGwgPSBodWxsLmZpbHRlcihmdW5jdGlvbihjZWxsKSB7XG4gICAgICB2YXIgY291bnQgPSAwXG4gICAgICBmb3IodmFyIGo9MDsgajw9ZDsgKytqKSB7XG4gICAgICAgIHZhciB2ID0gZGluZGV4W2NlbGxbal1dXG4gICAgICAgIGlmKHYgPCAwKSB7XG4gICAgICAgICAgaWYoKytjb3VudCA+PSAyKSB7XG4gICAgICAgICAgICByZXR1cm4gZmFsc2VcbiAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICAgICAgY2VsbFtqXSA9IHZcbiAgICAgIH1cbiAgICAgIHJldHVybiB0cnVlXG4gICAgfSlcbiAgfSBlbHNlIHtcbiAgICBodWxsID0gaHVsbC5maWx0ZXIoZnVuY3Rpb24oY2VsbCkge1xuICAgICAgZm9yKHZhciBpPTA7IGk8PWQ7ICsraSkge1xuICAgICAgICB2YXIgdiA9IGRpbmRleFtjZWxsW2ldXVxuICAgICAgICBpZih2IDwgMCkge1xuICAgICAgICAgIHJldHVybiBmYWxzZVxuICAgICAgICB9XG4gICAgICAgIGNlbGxbaV0gPSB2XG4gICAgICB9XG4gICAgICByZXR1cm4gdHJ1ZVxuICAgIH0pXG4gIH1cblxuICBpZihkICYgMSkge1xuICAgIGZvcih2YXIgaT0wOyBpPGh1bGwubGVuZ3RoOyArK2kpIHtcbiAgICAgIHZhciBoID0gaHVsbFtpXVxuICAgICAgdmFyIHggPSBoWzBdXG4gICAgICBoWzBdID0gaFsxXVxuICAgICAgaFsxXSA9IHhcbiAgICB9XG4gIH1cblxuICByZXR1cm4gaHVsbFxufSIsIlxubW9kdWxlLmV4cG9ydHMgPSBwYXJzZVxuXG4vKipcbiAqIGV4cGVjdGVkIGFyZ3VtZW50IGxlbmd0aHNcbiAqIEB0eXBlIHtPYmplY3R9XG4gKi9cblxudmFyIGxlbmd0aCA9IHthOiA3LCBjOiA2LCBoOiAxLCBsOiAyLCBtOiAyLCBxOiA0LCBzOiA0LCB0OiAyLCB2OiAxLCB6OiAwfVxuXG4vKipcbiAqIHNlZ21lbnQgcGF0dGVyblxuICogQHR5cGUge1JlZ0V4cH1cbiAqL1xuXG52YXIgc2VnbWVudCA9IC8oW2FzdHZ6cW1obGNdKShbXmFzdHZ6cW1obGNdKikvaWdcblxuLyoqXG4gKiBwYXJzZSBhbiBzdmcgcGF0aCBkYXRhIHN0cmluZy4gR2VuZXJhdGVzIGFuIEFycmF5XG4gKiBvZiBjb21tYW5kcyB3aGVyZSBlYWNoIGNvbW1hbmQgaXMgYW4gQXJyYXkgb2YgdGhlXG4gKiBmb3JtIGBbY29tbWFuZCwgYXJnMSwgYXJnMiwgLi4uXWBcbiAqXG4gKiBAcGFyYW0ge1N0cmluZ30gcGF0aFxuICogQHJldHVybiB7QXJyYXl9XG4gKi9cblxuZnVuY3Rpb24gcGFyc2UocGF0aCkge1xuXHR2YXIgZGF0YSA9IFtdXG5cdHBhdGgucmVwbGFjZShzZWdtZW50LCBmdW5jdGlvbihfLCBjb21tYW5kLCBhcmdzKXtcblx0XHR2YXIgdHlwZSA9IGNvbW1hbmQudG9Mb3dlckNhc2UoKVxuXHRcdGFyZ3MgPSBwYXJzZVZhbHVlcyhhcmdzKVxuXG5cdFx0Ly8gb3ZlcmxvYWRlZCBtb3ZlVG9cblx0XHRpZiAodHlwZSA9PSAnbScgJiYgYXJncy5sZW5ndGggPiAyKSB7XG5cdFx0XHRkYXRhLnB1c2goW2NvbW1hbmRdLmNvbmNhdChhcmdzLnNwbGljZSgwLCAyKSkpXG5cdFx0XHR0eXBlID0gJ2wnXG5cdFx0XHRjb21tYW5kID0gY29tbWFuZCA9PSAnbScgPyAnbCcgOiAnTCdcblx0XHR9XG5cblx0XHR3aGlsZSAodHJ1ZSkge1xuXHRcdFx0aWYgKGFyZ3MubGVuZ3RoID09IGxlbmd0aFt0eXBlXSkge1xuXHRcdFx0XHRhcmdzLnVuc2hpZnQoY29tbWFuZClcblx0XHRcdFx0cmV0dXJuIGRhdGEucHVzaChhcmdzKVxuXHRcdFx0fVxuXHRcdFx0aWYgKGFyZ3MubGVuZ3RoIDwgbGVuZ3RoW3R5cGVdKSB0aHJvdyBuZXcgRXJyb3IoJ21hbGZvcm1lZCBwYXRoIGRhdGEnKVxuXHRcdFx0ZGF0YS5wdXNoKFtjb21tYW5kXS5jb25jYXQoYXJncy5zcGxpY2UoMCwgbGVuZ3RoW3R5cGVdKSkpXG5cdFx0fVxuXHR9KVxuXHRyZXR1cm4gZGF0YVxufVxuXG5mdW5jdGlvbiBwYXJzZVZhbHVlcyhhcmdzKXtcblx0YXJncyA9IGFyZ3MubWF0Y2goLy0/Wy4wLTldKyg/OmVbLStdP1xcZCspPy9pZylcblx0cmV0dXJuIGFyZ3MgPyBhcmdzLm1hcChOdW1iZXIpIDogW11cbn1cbiIsIid1c2Ugc3RyaWN0JztcblxudmFyIHNpZ24gPSByZXF1aXJlKCcuL3V0aWxpdGllcy5qcycpLnNpZ247XG52YXIgY2FsY3VsYXRlRGlzdGFuY2UgPSByZXF1aXJlKCcuL3V0aWxpdGllcy5qcycpLmRpc3RhbmNlO1xuXG4vLyB2YXIgcG9pbnRzID0gcmVxdWlyZSgnLi9pbml0aWFsaXplUG9pbnRzLmpzJykucG9pbnRzO1xuLy8gdmFyIGNpdHlTZXQgPSByZXF1aXJlKCcuL2luaXRpYWxpemVQb2ludHMuanMnKS5jaXR5U2V0O1xuLy8gdmFyIHRleHRQb2ludHNJZCA9IHJlcXVpcmUoJy4vaW5pdGlhbGl6ZVBvaW50cy5qcycpLnRleHRQb2ludHNJZDtcbi8vIHZhciBwb3NzaWJsZVN0YXJ0UG9pbnRzSWQgPSByZXF1aXJlKCcuL2luaXRpYWxpemVQb2ludHMuanMnKS5wb3NzaWJsZVN0YXJ0UG9pbnRzSWQ7XG5cbnZhciBsaXZlTW91c2VQb3NpdGlvbiA9IHJlcXVpcmUoJy4vbW91c2UuanMnKTtcblxudmFyIFZlY3RvciA9IHJlcXVpcmUoJy4vdmVjdG9yLmpzJyk7XG5cbnZhciByYW5kb20gPSBNYXRoLnJhbmRvbTtcbnZhciBmbG9vciA9IE1hdGguZmxvb3I7XG4vLyB2YXIgUkVQVUxTSU9OID0gMC4wNTtcbi8vIHZhciBSRVBVTFNJT05TUEVFRCA9IDAuMDAyO1xuLy8gdmFyIEFOVFZFTE9DSVRZID0gMC4wMDE7XG5cbm1vZHVsZS5leHBvcnRzID0gZnVuY3Rpb24oY29udGFpbmVyLCBpbml0UG9pbnRzLCBvcHRpb25zKXtcblxuICAgIGNvbnNvbGUubG9nKCdPcHRpb25zIGFudCA6Jywgb3B0aW9ucyk7XG4gICAgLy8gRGVmaW5lIHRob3NlIHBhcmFtZXRlcnMgYXMgYXR0cmlidXRlcyBvZiBBbnQgb2JqZWN0ID9cbiAgICB2YXIgUkVQVUxTSU9OID0gb3B0aW9ucy5yZXBTaXplO1xuICAgIHZhciBSRVBVTFNJT05TUEVFRCA9IG9wdGlvbnMucmVwU3BlZWQ7XG4gICAgdmFyIEFOVFZFTE9DSVRZID0gb3B0aW9ucy52ZWxvY2l0eTtcbiAgICB2YXIgV0VJR0hUID0gb3B0aW9ucy53ZWlnaHQ7XG5cbiAgICB2YXIgbW91c2UgPSBsaXZlTW91c2VQb3NpdGlvbihjb250YWluZXIpO1xuXG4gICAgdmFyIHBvaW50cyA9IGluaXRQb2ludHMucG9pbnRzO1xuICAgIHZhciBjaXR5U2V0ID0gaW5pdFBvaW50cy5jaXR5U2V0O1xuICAgIHZhciB0ZXh0UG9pbnRzSWQgPSBpbml0UG9pbnRzLnRleHRQb2ludHNJZDtcbiAgICB2YXIgcG9zc2libGVTdGFydFBvaW50c0lkID0gaW5pdFBvaW50cy5wb3NzaWJsZVN0YXJ0UG9pbnRzSWQ7XG5cblxuICAgIGZ1bmN0aW9uIEFudChwb2ludCkge1xuICAgICAgICB0aGlzLnggPSBwb2ludC54OyAgICAgICAgICAgICAgICBcbiAgICAgICAgdGhpcy55ID0gcG9pbnQueTtcbiAgICAgICAgdGhpcy52ZWxvY2l0eSA9IEFOVFZFTE9DSVRZO1xuICAgICAgICB0aGlzLndlaWdodCA9IFdFSUdIVDtcbiAgICAgICAgdGhpcy5yZXBTaXplID0gUkVQVUxTSU9OO1xuICAgICAgICB0aGlzLnJlcFNwZWVkID0gUkVQVUxTSU9OU1BFRUQ7XG4gICAgICAgIHRoaXMuZWRnZSA9IHVuZGVmaW5lZDtcbiAgICAgICAgdGhpcy5zdGF0ZSA9IFwiZm9yYWdlXCI7XG4gICAgICAgIHRoaXMuZWRnZXMgPSBbXTtcbiAgICAgICAgdGhpcy5sYXN0Q2l0eSA9IHVuZGVmaW5lZDtcbiAgICAgICAgdGhpcy5vcmlnaW4gPSBwb2ludDtcbiAgICAgICAgdGhpcy5kZXN0aW5hdGlvbiA9IHVuZGVmaW5lZDtcbiAgICAgICAgdGhpcy5vcmllbnRhdGlvbiA9IHVuZGVmaW5lZDtcbiAgICAgICAgdGhpcy5kaXJlY3Rpb24gPSBuZXcgVmVjdG9yKDAsMCk7XG4gICAgICAgIHRoaXMucHJvZyA9IDA7XG4gICAgfVxuICAgIC8vIGZvcmFnZTogdGhlIGFudCB3YW5kZXJzIGFyb3VuZCB3aXRob3V0IGFueSBwaGVyb21vbiBkZXBvc2l0aW9uXG4gICAgLy8gb25jZSBpdCBmaW5kcyBhIGNpdHksIGl0IHN0YXJ0cyByZW1lbWJlcmluZyB0aGUgbm9kZXMgaXQgZ29lcyB0aHJvdWdoXG4gICAgLy8gd2hlbiBpdCBmaW5kcyBhbm90aGVyIGNpdHksIGl0IGNvbXB1dGVzIHRoZSBwYXRoIGxlbmd0aCBhbmQgYWRkcyBwaGVyb21vbnMgb25lIGVhY2ggZWRnZXNcbiAgICAvLyBwcm9wb3J0aW9ubmFseSB0byB0aGUgc2hvcnRlc3RuZXNzIG9mIHRoZSBwYXRoXG4gICAgLy8gaXQgcmVzZXRzIHRoZSBsaXN0IG9mIG5vZGVzIGFuZCBjb250aW51ZXNcbiAgICAvLyB3aGlsZSBmb3JhZ2luZyB0aGUgYW50IGNob3NlcyB0aGUgcGF0aCB3aXRoIGEgcGhlcm9tb24gcHJlZmVyZW5jZVxuXG5cbiAgICAvLyBzdGF0aWMgbWV0aG9kc1xuICAgIEFudC5nZW5lcmF0ZVJhbmRTdGFydFBvaW50ID0gZnVuY3Rpb24oKSB7XG4gICAgICAgIHZhciByYW5kSWQgPSBNYXRoLmZsb29yKHBvc3NpYmxlU3RhcnRQb2ludHNJZC5sZW5ndGggKiByYW5kb20oKSk7XG4gICAgICAgIHZhciByYW5kU3RhcnRQb2ludCA9IHBvaW50c1twb3NzaWJsZVN0YXJ0UG9pbnRzSWRbcmFuZElkXV07XG4gICAgICAgIHJldHVybiByYW5kU3RhcnRQb2ludDtcbiAgICB9XG5cblxuICAgIC8vIG1ldGhvZHNcbiAgICBBbnQucHJvdG90eXBlID0ge1xuXG4gICAgICAgIHRyYW5zaXQ6IGZ1bmN0aW9uKCl7XG4gICAgICAgICAgICBzd2l0Y2ggKHRoaXMuc3RhdGUpIHtcbiAgICAgICAgICAgIGNhc2UgXCJmb3JhZ2VcIjpcbiAgICAgICAgICAgICAgICB2YXIgcmVzID0gdGhpcy5tb3ZlKCk7XG4gICAgICAgICAgICAgICAgaWYgKHJlcy5jaXR5UmVhY2hlZCkge1xuICAgICAgICAgICAgICAgICAgICB0aGlzLnN0YXRlID0gXCJwaGVyb21vbmluZ1wiO1xuICAgICAgICAgICAgICAgICAgICB0aGlzLmxhc3RDaXR5ID0gdGhpcy5vcmlnaW4uaWQ7XG4gICAgICAgICAgICAgICAgfTtcbiAgICAgICAgICAgICAgICBicmVhaztcbiAgICAgICAgICAgIGNhc2UgXCJwaGVyb21vbmluZ1wiOlxuICAgICAgICAgICAgICAgIHZhciByZXMgPSB0aGlzLm1vdmUoKTtcbiAgICAgICAgICAgICAgICBpZiAocmVzLmVkZ2VDaGFuZ2VkKSB7XG4gICAgICAgICAgICAgICAgICAgIHRoaXMuZWRnZXMucHVzaCh0aGlzLmVkZ2UpO1xuICAgICAgICAgICAgICAgICAgICAvLyBmb3VuZCBhIGNpdHlcbiAgICAgICAgICAgICAgICAgICAgaWYgKHJlcy5jaXR5UmVhY2hlZCAmJiAodGhpcy5vcmlnaW4uaWQgIT0gdGhpcy5sYXN0Q2l0eSkgKXtcbiAgICAgICAgICAgICAgICAgICAgICAgIC8vIGNvbXB1dGUgdGhlIGxlbmd0aCBvZiB0aGUgcGF0aFxuICAgICAgICAgICAgICAgICAgICAgICAgdmFyIHBhdGhMZW5ndGggPSB0aGlzLmVkZ2VzLm1hcChmdW5jdGlvbihlKXtyZXR1cm4gZS5kaXN0YW5jZX0pLnJlZHVjZShmdW5jdGlvbihhLGIpe3JldHVybiBhICsgYn0pO1xuICAgICAgICAgICAgICAgICAgICAgICAgdmFyIGRlbHRhUGhlcm9tb25lID0gMS9wYXRoTGVuZ3RoO1xuICAgICAgICAgICAgICAgICAgICAgICAgdmFyIGFudFdlaWdodCA9IHRoaXMud2VpZ2h0O1xuICAgICAgICAgICAgICAgICAgICAgICAgdGhpcy5lZGdlcy5mb3JFYWNoKGZ1bmN0aW9uKGUpe1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgIHZhciBhID0gZS5wdDEsIGIgPSBlLnB0Miwgd2VpZ2h0ID0gMTsgIFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgIC8vIGluY3JlYXNlZCBkcm9wcGVkIHBoZXJvbW9ucyBmb3IgdGV4dEVkZ2VzXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgaWYgKChjaXR5U2V0LmluZGV4T2YoYS5pZCkgIT0gLTEpICYmIGNpdHlTZXQuaW5kZXhPZihiLmlkKSAhPSAtMSAmJiAoTWF0aC5hYnMoYS5pZCAtIGIuaWQpID09IDEpKVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgd2VpZ2h0ICo9IGFudFdlaWdodDtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgZS5waGVyb21vbiArPSAoZGVsdGFQaGVyb21vbmUgKiB3ZWlnaHQpO1xuICAgICAgICAgICAgICAgICAgICAgICAgfSk7XG5cbiAgICAgICAgICAgICAgICAgICAgICAgIHRoaXMuZWRnZXMgPSBbdGhpcy5lZGdlXTtcbiAgICAgICAgICAgICAgICAgICAgICAgIHRoaXMubGFzdENpdHkgPSB0aGlzLm9yaWdpbi5pZDtcbiAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgYnJlYWs7XG4gICAgICAgICAgICB9XG5cbiAgICAgICAgfSxcblxuICAgICAgICBzZXREaXJlY3Rpb246IGZ1bmN0aW9uKCl7XG4gICAgICAgICAgICB2YXIgcG9zc2libGVFZGdlcyA9IFtdO1xuXG4gICAgICAgICAgICBmb3IgKHZhciBpID0gMDsgaSA8IHRoaXMub3JpZ2luLm5leHRzLmxlbmd0aDsgaSsrKVxuICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgIHBvc3NpYmxlRWRnZXNbaV0gPSB0aGlzLm9yaWdpbi5uZXh0c1tpXTtcbiAgICAgICAgICAgIH0gXG5cbiAgICAgICAgICAgIC8vIGNvbnNvbGUubG9nKCdzbWVsbHMxOiAnLCBwb3NzaWJsZUVkZ2VzKTtcblxuICAgICAgICAgICAgcG9zc2libGVFZGdlcy5zcGxpY2UocG9zc2libGVFZGdlcy5pbmRleE9mKHRoaXMuZWRnZSksMSk7XG5cbiAgICAgICAgICAgIC8vIGZsaXAgYSBjb2luIGFuZCBlaXRoZXIgdGFrZSB0aGUgc21lbGxpZXN0IHBhdGggb3IgYSByYW5kb20gb25lXG4gICAgICAgICAgICBpZiAocmFuZG9tKCkgPiAwLjUpe1xuICAgICAgICAgICAgICAgIHZhciBzbWVsbHMgPSBwb3NzaWJsZUVkZ2VzLm1hcChmdW5jdGlvbihlKXtyZXR1cm4gZS5waGVyb21vbjt9KTtcbiAgICAgICAgICAgICAgICB2YXIgaW5kZXggPSBzbWVsbHMuaW5kZXhPZihNYXRoLm1heC5hcHBseShNYXRoLCBzbWVsbHMpKTtcbiAgICAgICAgICAgICAgICB0aGlzLmVkZ2UgPSBwb3NzaWJsZUVkZ2VzW2luZGV4XTtcbiAgICAgICAgICAgIH0gXG4gICAgICAgICAgICBlbHNle1xuICAgICAgICAgICAgICAgIHRoaXMuZWRnZSA9IHBvc3NpYmxlRWRnZXNbZmxvb3IocmFuZG9tKCkqcG9zc2libGVFZGdlcy5sZW5ndGgpXTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICBcblxuICAgICAgICAgICAgLy8gc2V0IHRoZSBkZXN0aW5hdGlvbiBwb2ludCwgYmVpbmcgZWRnZS5wdDEgb3IgZWRnZS5wdDJcbiAgICAgICAgICAgIHRoaXMuZGVzdGluYXRpb24gPSAodGhpcy5vcmlnaW4gPT0gdGhpcy5lZGdlLnB0MSkgPyB0aGlzLmVkZ2UucHQyIDogdGhpcy5lZGdlLnB0MTtcblxuICAgICAgICAgICAgdGhpcy5kaXJlY3Rpb24ueCA9IHRoaXMuZGVzdGluYXRpb24ueCAtIHRoaXMub3JpZ2luLng7IFxuICAgICAgICAgICAgdGhpcy5kaXJlY3Rpb24ueSA9IHRoaXMuZGVzdGluYXRpb24ueSAtIHRoaXMub3JpZ2luLnk7XG5cbiAgICAgICAgICAgIHRoaXMuZGlyZWN0aW9uLm5vcm1hbGl6ZSgpO1xuICAgICAgICB9LFxuXG4gICAgICAgIG1vdmU6IGZ1bmN0aW9uKCl7XG4gICAgICAgICAgICAvLyBjb25zb2xlLmxvZygnbW92ZScpO1xuICAgICAgICAgICAgdmFyIGVkZ2VDaGFuZ2VkO1xuICAgICAgICAgICAgdmFyIGNpdHlSZWFjaGVkID0gZmFsc2U7XG5cbiAgICAgICAgICAgIHRoaXMuZGlyZWN0aW9uLnggPSB0aGlzLmRlc3RpbmF0aW9uLnggLSB0aGlzLng7IFxuICAgICAgICAgICAgdGhpcy5kaXJlY3Rpb24ueSA9IHRoaXMuZGVzdGluYXRpb24ueSAtIHRoaXMueTtcbiAgICAgICAgICAgIHRoaXMuZGlyZWN0aW9uLm5vcm1hbGl6ZSgpO1xuXG4gICAgICAgICAgICAvLyBvbiBlZGdlXG4gICAgICAgICAgICBpZiAoKGNhbGN1bGF0ZURpc3RhbmNlKHRoaXMsIHRoaXMuZGVzdGluYXRpb24pID4gdGhpcy5yZXBTcGVlZCkpe1xuXG4gICAgICAgICAgICAgICAgLy8gYSBkZWx0YSBtb3ZlbWVudCB3aWxsIGJlIGFwcGxpZWQgaWYgY29sbGlzaW9uIHdpdGggb2JzdGFjbGUgZGV0ZWN0ZWRcbiAgICAgICAgICAgICAgICB2YXIgZGVsdGEgPSB0aGlzLmF2b2lkT2JzdGFjbGUoKTtcblxuICAgICAgICAgICAgICAgIHRoaXMueCArPSB0aGlzLnZlbG9jaXR5ICogdGhpcy5kaXJlY3Rpb24ueCArIGRlbHRhLnggKiB0aGlzLnJlcFNwZWVkO1xuICAgICAgICAgICAgICAgIHRoaXMueSArPSB0aGlzLnZlbG9jaXR5ICogdGhpcy5kaXJlY3Rpb24ueSArIGRlbHRhLnkgKiB0aGlzLnJlcFNwZWVkO1xuXG4gICAgICAgICAgICAgICAgdGhpcy5wcm9nID0gdGhpcy5jYWxjdWxhdGVQcm9ncmVzc2lvbigpO1xuICAgICAgICAgICAgICAgIFxuICAgICAgICAgICAgICAgIGVkZ2VDaGFuZ2VkID0gZmFsc2U7XG5cbiAgICAgICAgICAgIC8vIG9uIHZlcnRleFxuICAgICAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgICAgICAvLyBjb25zb2xlLmxvZygncmVhY2hlZCcpO1xuICAgICAgICAgICAgICAgIHRoaXMuc3RlcCA9IDA7XG4gICAgICAgICAgICAgICAgdGhpcy5wcm9nID0gMDtcbiAgICAgICAgICAgICAgICB0aGlzLm9yaWdpbiA9IHRoaXMuZGVzdGluYXRpb247XG4gICAgICAgICAgICAgICAgdGhpcy54ID0gdGhpcy5vcmlnaW4ueDtcbiAgICAgICAgICAgICAgICB0aGlzLnkgPSB0aGlzLm9yaWdpbi55O1xuXG4gICAgICAgICAgICAgICAgdGhpcy5zZXREaXJlY3Rpb24oKTtcblxuICAgICAgICAgICAgICAgIGNpdHlSZWFjaGVkID0gKGNpdHlTZXQuaW5kZXhPZih0aGlzLm9yaWdpbi5pZCkgIT0gLTEpO1xuICAgICAgICAgICAgICAgIGVkZ2VDaGFuZ2VkID0gdHJ1ZTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgICAgIHJldHVybiB7Y2l0eVJlYWNoZWQ6IGNpdHlSZWFjaGVkLCBlZGdlQ2hhbmdlZDogZWRnZUNoYW5nZWR9O1xuICAgICAgICB9LFxuXG4gICAgICAgIGF2b2lkT2JzdGFjbGU6IGZ1bmN0aW9uKCl7XG4gICAgICAgICAgICB2YXIgZGlzdGFuY2UgPSBjYWxjdWxhdGVEaXN0YW5jZSh0aGlzLCBtb3VzZSk7XG4gICAgICAgIFxuICAgICAgICAgICAgaWYgKGRpc3RhbmNlIDw9IHRoaXMucmVwU2l6ZSkge1xuXG4gICAgICAgICAgICAgICAgcmV0dXJuIHtcbiAgICAgICAgICAgICAgICAgICAgLy8gZGVsdGEgbW92ZW1lbnQgaXMgY29tcG9zZWQgb2YgYSByZXB1bHNpb24gZGVsdGEgYW5kIGEgY2lyY3VsYXIgZGVsdGEgXG4gICAgICAgICAgICAgICAgICAgIHg6ICh0aGlzLnggLSBtb3VzZS54KS9kaXN0YW5jZSArICh0aGlzLnkgLSBtb3VzZS55KS9kaXN0YW5jZSAqIDEsXG4gICAgICAgICAgICAgICAgICAgIHk6ICh0aGlzLnkgLSBtb3VzZS55KS9kaXN0YW5jZSAtICh0aGlzLnggLSBtb3VzZS54KS9kaXN0YW5jZSAqIDFcbiAgICAgICAgICAgICAgICB9O1xuICAgICAgICAgICAgfVxuICAgICAgICAgICAgZWxzZVxuICAgICAgICAgICAgICAgIHJldHVybiB7eDowLCB5OjB9O1xuICAgICAgICB9LFxuXG4gICAgICAgIGNhbGN1bGF0ZVByb2dyZXNzaW9uOiBmdW5jdGlvbigpe1xuICAgICAgICAgICAgdmFyIHYgPSBuZXcgVmVjdG9yKHRoaXMueCAtIHRoaXMub3JpZ2luLngsIHRoaXMueSAtIHRoaXMub3JpZ2luLnkpO1xuICAgICAgICAgICAgdmFyIG5vcm0gPSB2Lm5vcm0oKTtcblxuICAgICAgICAgICAgdmFyIHRoZXRhID0gKHYueCAqIHRoaXMuZWRnZS5kaXJlY3Rpb24ueCArIHYueSAqIHRoaXMuZWRnZS5kaXJlY3Rpb24ueSkgLyBub3JtO1xuICAgICAgICAgICAgdmFyIHByb2cgPSBub3JtICogTWF0aC5hYnModGhldGEpO1xuICAgICAgICAgICAgLy8gcmV0dXJucyBsZW5ndGggb2YgcHJvamVjdGlvbiBvbiBlZGdlXG4gICAgICAgICAgICByZXR1cm4gcHJvZztcbiAgICAgICAgfVxuXG4gICAgfTtcbiAgICByZXR1cm4gQW50O1xufVxuXG4iLCIndXNlIHN0cmljdCdcblxubW9kdWxlLmV4cG9ydHMgPSBmdW5jdGlvbiAoQW50KSB7XG5cblx0dmFyIG5iQW50c1BlclN0ZXAgPSAxMDA7XG5cblx0ZnVuY3Rpb24gY3JlYXRlR3JvdXAocG9wdWxhdGlvbil7XG5cdFx0Zm9yICh2YXIgaSA9IDA7IGkgPCBuYkFudHNQZXJTdGVwOyBpKyspIHtcblx0XHRcdHZhciBuZXdBbnQgPSBuZXcgQW50KEFudC5nZW5lcmF0ZVJhbmRTdGFydFBvaW50KCkpO1xuXHRcdFx0bmV3QW50LnNldERpcmVjdGlvbigpO1xuXHRcdFx0cG9wdWxhdGlvbi5wdXNoKG5ld0FudCk7XG5cdFx0fVxuXG4vLyBcdFx0Y29uc29sZS5sb2coJ0NyZWF0ZWQgQW50cyBHcm91cDogXFxcbi8vICgrICcgKyBuYkFudHNQZXJTdGVwICsgJykgPT4gJyArIHBvcHVsYXRpb24ubGVuZ3RoKTtcblxuXHRcdHJldHVybiBwb3B1bGF0aW9uO1xuXHR9XG5cblx0ZnVuY3Rpb24gcmVtb3ZlR3JvdXAocG9wdWxhdGlvbiwgbmJEZWFkKXtcblx0XHRwb3B1bGF0aW9uID0gcG9wdWxhdGlvbi5zbGljZSgwLCBwb3B1bGF0aW9uLmxlbmd0aCAtIG5iRGVhZCk7XG5cbi8vIFx0XHRjb25zb2xlLmxvZygnUmVtb3ZlZCBBbnRzIEdyb3VwOiBcXFxuLy8gKC0gJyArIG5iQW50c1BlclN0ZXAgKyAnKSA9PiAnICsgcG9wdWxhdGlvbi5sZW5ndGgpO1xuXG5cdFx0cmV0dXJuIHBvcHVsYXRpb247XG5cblx0fVxuXG5cdHJldHVybiB7XG5cdFx0Y3JlYXRlOiBjcmVhdGVHcm91cCxcblx0XHRyZW1vdmU6IHJlbW92ZUdyb3VwXG5cdH07XG5cbn1cblx0IiwiJ3VzZSBzdHJpY3QnXG5cbnZhciBkdCA9IHJlcXVpcmUoXCJkZWxhdW5heS10cmlhbmd1bGF0ZVwiKTtcblxudmFyIHJhbmdlID0gcmVxdWlyZSgnLi91dGlsaXRpZXMuanMnKS5yYW5nZTtcblxuLy8gdmFyIHBvaW50cyA9IHJlcXVpcmUoJy4vaW5pdGlhbGl6ZVBvaW50cy5qcycpLnBvaW50cztcbnZhciB0ZXh0TWVzaCA9IHJlcXVpcmUoJy4vaW5pdGlhbGl6ZVBvaW50cy5qcycpLnRleHRNZXNoO1xudmFyIGNpdHlTZXQgPSByZXF1aXJlKCcuL2luaXRpYWxpemVQb2ludHMuanMnKS5jaXR5U2V0O1xudmFyIG5iUmFuZG9tUG9pbnRzID0gcmVxdWlyZSgnLi9pbml0aWFsaXplUG9pbnRzLmpzJykubmJSYW5kb21Qb2ludHM7XG52YXIgZm9yY2VkRWRnZXMgPSByZXF1aXJlKCcuL2luaXRpYWxpemVQb2ludHMuanMnKS5mb3JjZWRFZGdlcztcblxudmFyIEVkZ2UgPSByZXF1aXJlKCcuL2VkZ2UuanMnKTtcblxuXG5tb2R1bGUuZXhwb3J0cyA9IGZ1bmN0aW9uKHBvaW50cyl7XG4gICAgLy8gdHJpYW5ndWxhdGVcbiAgICB2YXIgY2VsbHMgPSBkdChwb2ludHMubWFwKGZ1bmN0aW9uKHApe1xuICAgICAgICByZXR1cm4gW3AueCwgcC55XTtcbiAgICB9KSk7XG5cbiAgICB2YXIgZWRnZXMgPSBbXTtcbiAgICB2YXIgcGVybXV0YXRpb25zID0gW1swLDFdLCBbMCwyXSwgWzEsMl1dO1xuXG4gICAgLy8gZm9yY2UgdGhlIGVkZ2VzIG9mIHRoZSB0ZXh0IHRvIGJlIGVkZ2VzIG9mIHRoZSBncmFwaFxuICAgIGlmICh0ZXh0TWVzaCkge1xuICAgICAgICByYW5nZSgwLCBwb2ludHMubGVuZ3RoIC0gbmJSYW5kb21Qb2ludHMpLmZvckVhY2goZnVuY3Rpb24oaWQpe1xuICAgICAgICAgICAgdmFyIGRpcmVjdExpbmsgPSBmb3JjZWRFZGdlc1tpZF07XG4gICAgICAgICAgICB2YXIgdGV4dEVkZ2UgPSBFZGdlLmNyZWF0ZShwb2ludHNbaWRdLCBwb2ludHNbZGlyZWN0TGlua10pO1xuICAgICAgICAgICAgZWRnZXMucHVzaCh0ZXh0RWRnZSk7XG4gICAgICAgICAgICBwb2ludHNbaWRdLm5leHRzLnB1c2godGV4dEVkZ2UpO1xuICAgICAgICB9KVxuICAgIH1cblxuXG4gICAgY2VsbHMuZm9yRWFjaChmdW5jdGlvbihjZWxsKXtcbiAgICAgICBcbiAgICAgICAgZm9yICh2YXIgaSA9IDA7IGkgPCAzOyArK2kpeyAgLy8gZm9yIGVhY2ggcG9pbnQuaWQgbGlzdGVkIGluIGN1cnJlbnQgY2VsbFxuICAgICAgICAgICAgdmFyIHB0ID0gcG9pbnRzW2NlbGxbaV1dO1xuXG4gICAgICAgICAgICBmb3IgKHZhciBqID0gMTsgaiA8IDM7ICsrail7IFxuXG4gICAgICAgICAgICAgICAgdmFyIHB0aiA9IHBvaW50c1tjZWxsWyggaSArIGogKSAlIDNdXTsgLy8gcGljayBvbmUgb2YgdGhlIG90aGVyIDIgcG9pbnRzIG9mIHRoZSBjZWxsXG4gICAgICAgICAgICAgICAgdmFyIG5ld0VkZ2UgPSB1bmRlZmluZWQ7XG5cbiAgICAgICAgICAgICAgICAvLyBpZiBwdCBhbHJlYWR5IGhhcyBuZXh0RWRnZXMgLi4uXG4gICAgICAgICAgICAgICAgaWYgKHB0Lm5leHRzLmxlbmd0aCAhPSAwKSB7XG4gICAgICAgICAgICAgICAgICAgIFxuICAgICAgICAgICAgICAgICAgICAvLyAuLi4gZ2V0IHRoZSBwb2ludHMgY29ycmVzcG9uZGluZyAuLi5cbiAgICAgICAgICAgICAgICAgICAgdmFyIHRlbXBQb2ludHMgPSBwdC5uZXh0cy5tYXAoZnVuY3Rpb24oZSl7XG4gICAgICAgICAgICAgICAgICAgICAgICByZXR1cm4gW2UucHQxLCBlLnB0Ml07XG4gICAgICAgICAgICAgICAgICAgIH0pLnJlZHVjZShmdW5jdGlvbihhLCBiKXtcbiAgICAgICAgICAgICAgICAgICAgICAgICByZXR1cm4gYS5jb25jYXQoYik7XG4gICAgICAgICAgICAgICAgICAgIH0pO1xuXG4gICAgICAgICAgICAgICAgICAgIC8vIC4uLiBhbmQgY2hlY2sgaWYgcHRqIGFscmVhZHkgaXMgcGFydCBvZiB0aGUgZXhpc3RpbmcgbmV4dEVkZ2VzLiBJZiBub3QsIGFkZCB0aGUgZWRnZS5cbiAgICAgICAgICAgICAgICAgICAgaWYgKHRlbXBQb2ludHMuaW5kZXhPZihwdGopID09IC0xKXtcbiAgICAgICAgICAgICAgICAgICAgICAgIG5ld0VkZ2UgPSBFZGdlLmNyZWF0ZShwdCwgcHRqKTtcbiAgICAgICAgICAgICAgICAgICAgICAgIGVkZ2VzLnB1c2gobmV3RWRnZSk7XG4gICAgICAgICAgICAgICAgICAgICAgICBwdC5uZXh0cy5wdXNoKG5ld0VkZ2UpO1xuICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIGVsc2Uge1xuICAgICAgICAgICAgICAgICAgICBuZXdFZGdlID0gRWRnZS5jcmVhdGUocHQsIHB0aik7XG4gICAgICAgICAgICAgICAgICAgIGVkZ2VzLnB1c2gobmV3RWRnZSk7XG4gICAgICAgICAgICAgICAgICAgIHB0Lm5leHRzLnB1c2gobmV3RWRnZSk7XG4gICAgICAgICAgICAgICAgfVxuXG4gICAgICAgICAgICAgICAgLy8gYWRkIGFsc28gdGhlIGVkZ2UgdG8gdGhlIGVkZ2UncyBvdGhlciBwb2ludCdzIG5leHRFZGdlc1xuICAgICAgICAgICAgICAgIGlmIChuZXdFZGdlICE9IHVuZGVmaW5lZCl7XG4gICAgICAgICAgICAgICAgICAgIHB0ai5uZXh0cy5wdXNoKG5ld0VkZ2UpO1xuICAgICAgICAgICAgICAgIH0gICAgICAgICBcbiAgICAgICAgICAgIH1cblxuICAgICAgICAgICAgLy8gYWRkIHRoZSB0ZXh0RWRnZXMgdG8gbmV4dEVkZ2VzIG1hcFxuICAgICAgICAgICAgaWYgKHRleHRNZXNoICYmIChjaXR5U2V0LmluZGV4T2YocHQpICE9IC0xKSkge1xuICAgICAgICAgICAgICAgIHZhciB0ZXh0RWRnZSA9IEVkZ2UuY3JlYXRlKHB0LCBwb2ludHNbcHQuaWQgKyAxXSk7XG4gICAgICAgICAgICAgICAgZWRnZXMucHVzaCh0ZXh0RWRnZSk7XG4gICAgICAgICAgICAgICAgcHQubmV4dHMucHVzaCh0ZXh0RWRnZSk7XG4gICAgICAgICAgICB9XG5cbiAgICAgICAgfVxuICAgIH0pO1xuXG4gICAgcmV0dXJuIGVkZ2VzO1xufTsiLCIndXNlIHN0cmljdCc7XG5cbnZhciBzcXJ0ID0gTWF0aC5zcXJ0O1xudmFyIHBvdyA9IE1hdGgucG93O1xudmFyIGFicyA9IE1hdGguYWJzO1xudmFyIGF0YW4gPSBNYXRoLmF0YW47XG5cbnZhciBWZWN0b3IgPSByZXF1aXJlKCcuL3ZlY3Rvci5qcycpO1xuXG5cbmZ1bmN0aW9uIEVkZ2UocHRBLCBwdEIpIHtcbiAgICB2YXIgZGlzdGFuY2UgPSBzcXJ0KCBwb3cocHRBLnggLSBwdEIueCwgMikgKyBwb3cocHRBLnkgLSBwdEIueSwgMikgKTtcblxuICAgIC8vIGZpbmQgbGluZSBlcXVhdGlvbiBheCArIGJ5ICsgYyA9IDBcbiAgICB2YXIgYSA9IDE7XG4gICAgdmFyIGIgPSAtIChwdEIueCAtIHB0QS54KSAvIChwdEIueSAtIHB0QS55KTtcblxuICAgIC8vIG9yaWVudGF0ZSB2ZWN0b3IgKGEsYilcbiAgICBpZiAoYiA8IDApe1xuICAgICAgICBiID0gLWI7XG4gICAgICAgIGEgPSAtYTtcbiAgICB9XG5cbiAgICAvLyBub3JtYWxpemUgdmVjdG9yIChhLGIpXG4gICAgdmFyIG4gPSBuZXcgVmVjdG9yKGEsIGIpO1xuICAgIG4ubm9ybWFsaXplKCk7XG5cbiAgICB2YXIgYyA9IC0gKGEgKiBwdEEueCArIGIgKiBwdEEueSk7XG5cbiAgICAvLyAvLyBjYWxjdWxhdGUgdmVjdG9yIGRpcmVjdG9yXG4gICAgdmFyIHYgPSBuZXcgVmVjdG9yKHB0Qi54IC0gcHRBLngsIHB0Qi55IC0gcHRBLnkpO1xuICAgIFxuICAgIHYubm9ybWFsaXplKCk7XG5cbiAgICB0aGlzLmlkID0gdW5kZWZpbmVkO1xuICAgIHRoaXMucHQxID0gcHRBO1xuICAgIHRoaXMucHQyID0gcHRCO1xuICAgIHRoaXMuZGlyZWN0aW9uID0gdjtcbiAgICB0aGlzLm9ydGhEaXJlY3Rpb24gPSBuOyBcbiAgICB0aGlzLmRpc3RhbmNlID0gZGlzdGFuY2U7XG4gICAgdGhpcy5waGVyb21vbiA9IDEvZGlzdGFuY2U7XG4gICAgdGhpcy5saW5lID0ge1xuICAgICAgICBhOiBhLFxuICAgICAgICBiOiBiLFxuICAgICAgICBjOiBjLFxuICAgIH07XG5cbiAgICBpZiAodGhpcy5kaXN0YW5jZSA9PT0gMCkgY29uc29sZS5sb2coJ1pFUk8gIScpO1xufVxuXG5cbi8vIHN0YXRpYyBtZXRob2RzXG5FZGdlLmNyZWF0ZSA9IGZ1bmN0aW9uKHB0QSwgcHRCKSB7XG4gICAgdmFyIGVkZ2UgPSBuZXcgRWRnZShwdEEsIHB0Qik7XG4gICAgcmV0dXJuIGVkZ2U7XG59XG5cblxuLy8gbWV0aG9kc1xuRWRnZS5wcm90b3R5cGUgPSB7XG5cbiAgICBnZXRPdGhlclBvaW50OiBmdW5jdGlvbihwb2ludCkge1xuICAgICAgICBpZiAocG9pbnQgPT0gdGhpcy5wdDEpXG4gICAgICAgICAgICByZXR1cm4gdGhpcy5wdDI7XG4gICAgICAgIGVsc2UgaWYgKHBvaW50ID09IHRoaXMucHQyKVxuICAgICAgICAgICAgcmV0dXJuIHRoaXMucHQxO1xuICAgICAgICBlbHNlXG4gICAgICAgICAgICBjb25zb2xlLmxvZyhcIkVycm9yXCIpO1xuICAgIH0sXG5cbiAgICBjYWxjdWxhdGVEaXN0YW5jZTogZnVuY3Rpb24oeCwgeSkge1xuICAgICAgICB2YXIgYSA9IHRoaXMubGluZS5hLFxuICAgICAgICAgICAgYiA9IHRoaXMubGluZS5iLFxuICAgICAgICAgICAgYyA9IHRoaXMubGluZS5jO1xuICAgICAgICByZXR1cm4gYWJzKGEgKiB4ICsgYiAqIHkgKyBjKSAvIE1hdGguc3FydChNYXRoLnBvdyhhLDIpICsgTWF0aC5wb3coYiwyKSk7XG4gICAgfSxcblxufVxubW9kdWxlLmV4cG9ydHMgPSBFZGdlOyIsIid1c2Ugc3RyaWN0JztcblxudmFyIHBhcnNlID0gcmVxdWlyZSgncGFyc2Utc3ZnLXBhdGgnKTtcblxudmFyIHJhbmdlID0gcmVxdWlyZSgnLi91dGlsaXRpZXMuanMnKS5yYW5nZTtcblxudmFyIFBvaW50ID0gcmVxdWlyZSgnLi9wb2ludC5qcycpO1xudmFyIHN2Z1BhdGggPSByZXF1aXJlKCcuL3N2Z1BhdGguanMnKTtcblxudmFyIHJhbmRvbSA9IE1hdGgucmFuZG9tO1xuXG52YXIgbmJDaXR5ID0gMjtcblxudmFyIHRleHRNZXNoID0gdHJ1ZTtcblxuLy8gRnJhbWUgZGVmaW5pdGlvblxudmFyIHhJbml0ID0gMCwgeUluaXQgPSAwO1xudmFyIHcgPSAxLFxuICAgIGggPSAxO1xuXG52YXIgc3ZnU3RyaW5nID0gc3ZnUGF0aDtcblxuZnVuY3Rpb24gc3ZnVG9Qb2ludHMoc3ZnU3RyaW5nKSB7XG4gICAgdmFyIHBvaW50cyA9IFtdO1xuICAgIHZhciBlZGdlcyA9IE9iamVjdC5jcmVhdGUobnVsbCk7XG5cbiAgICB2YXIgYmVnaW5pbmdQYXRoO1xuXG4gICAgdmFyIFggPSAwO1xuICAgIHZhciBZID0gMDtcbiAgICB2YXIgbmJQb2ludHMgPSAwO1xuICAgIHZhciBwcmV2UG9pbnQ7XG5cbiAgICB2YXIgY29tbWFuZHMgPSBwYXJzZShzdmdTdHJpbmcpXG4gICAgZm9yICh2YXIgaT0wOyBpPGNvbW1hbmRzLmxlbmd0aDsgaSsrKXtcbiAgICAgICAgdmFyIGNvbW1hbmQgPSBjb21tYW5kc1tpXTtcbiAgICAgICAgc3dpdGNoIChjb21tYW5kWzBdKSB7XG4gICAgICAgICAgICBjYXNlIFwibVwiOlxuICAgICAgICAgICAgICAgIFggKz0gY29tbWFuZFsxXTtcbiAgICAgICAgICAgICAgICBZICs9IGNvbW1hbmRbMl07XG4gICAgICAgICAgICAgICAgcHJldlBvaW50ID0gdW5kZWZpbmVkO1xuICAgICAgICAgICAgICAgIGJlZ2luaW5nUGF0aCA9IG5iUG9pbnRzO1xuICAgICAgICAgICAgICAgIGJyZWFrO1xuICAgICAgICAgICAgY2FzZSBcIk1cIjpcbiAgICAgICAgICAgICAgICBYID0gY29tbWFuZFsxXTtcbiAgICAgICAgICAgICAgICBZID0gY29tbWFuZFsyXTtcbiAgICAgICAgICAgICAgICBwcmV2UG9pbnQgPSB1bmRlZmluZWQ7XG4gICAgICAgICAgICAgICAgYmVnaW5pbmdQYXRoID0gbmJQb2ludHM7XG4gICAgICAgICAgICAgICAgYnJlYWs7IFxuICAgICAgICAgICAgY2FzZSBcImNcIjpcbiAgICAgICAgICAgICAgICBYID0gY29tbWFuZFsxXTtcbiAgICAgICAgICAgICAgICBZID0gY29tbWFuZFsyXTtcbiAgICAgICAgICAgICAgICBwb2ludHMucHVzaCh7aWQ6bmJQb2ludHMsIHg6WCwgeTpZfSk7XG5cbiAgICAgICAgICAgICAgICBpZiAocHJldlBvaW50ICE9IHVuZGVmaW5lZCkge1xuICAgICAgICAgICAgICAgICAgICBlZGdlc1twcmV2UG9pbnRdID0gbmJQb2ludHM7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIHByZXZQb2ludCA9IG5iUG9pbnRzO1xuICAgICAgICAgICAgICAgIG5iUG9pbnRzKys7XG4gICAgICAgICAgICAgICAgYnJlYWs7IFxuICAgICAgICAgICAgY2FzZSBcImxcIjpcbiAgICAgICAgICAgICAgICBYID0gY29tbWFuZFsxXTtcbiAgICAgICAgICAgICAgICBZID0gY29tbWFuZFsyXTtcbiAgICAgICAgICAgICAgICBwb2ludHMucHVzaCh7aWQ6bmJQb2ludHMsIHg6WCwgeTpZfSk7XG5cbiAgICAgICAgICAgICAgICBpZiAocHJldlBvaW50ICE9IHVuZGVmaW5lZCkge1xuICAgICAgICAgICAgICAgICAgICBlZGdlc1twcmV2UG9pbnRdID0gbmJQb2ludHM7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIHByZXZQb2ludCA9IG5iUG9pbnRzO1xuICAgICAgICAgICAgICAgIG5iUG9pbnRzKys7XG4gICAgICAgICAgICAgICAgYnJlYWs7XG4gICAgICAgICAgICBjYXNlIFwielwiOlxuICAgICAgICAgICAgICAgIGVkZ2VzW3ByZXZQb2ludF0gPSBuYlBvaW50cztcbiAgICAgICAgICAgICAgICBiZWdpbmluZ1BhdGggPSB1bmRlZmluZWQ7XG4gICAgICAgICAgICAgICAgcHJldlBvaW50ID0gdW5kZWZpbmVkO1xuICAgICAgICAgICAgICAgIGJyZWFrOyAgICBcbiAgICAgICAgfVxuICAgIH1cbiAgICByZXR1cm4ge3BvaW50cyA6IHBvaW50cywgZWRnZXMgOiBlZGdlc307XG59XG5cbi8vIGluaXRpYWxpemUgcG9pbnRzXG5cbm1vZHVsZS5leHBvcnRzID0gZnVuY3Rpb24obmJTdGFydFBvaW50cywgbmJSYW5kb21Qb2ludHMpe1xuICAgIHZhciBwb2ludHMgPSBbXTtcbiAgICB2YXIgZm9yY2VkRWRnZXM7XG4gICAgdmFyIGNpdHlTZXQ7XG5cbiAgICBpZiAodGV4dE1lc2gpe1xuXG4gICAgICAgIHZhciBteVRleHQgPSBzdmdUb1BvaW50cyhzdmdTdHJpbmcpO1xuICAgICAgICBwb2ludHMgPSBteVRleHQucG9pbnRzO1xuICAgICAgICBmb3JjZWRFZGdlcyA9IG15VGV4dC5lZGdlcztcbiAgICAgICAgY2l0eVNldCA9IHJhbmdlKDAsIHBvaW50cy5sZW5ndGgpO1xuXG4gICAgICAgIGNvbnNvbGUubG9nKHBvaW50cyk7XG5cbiAgICAgICAgdmFyIHNjYWxlWCA9IDAuNTtcbiAgICAgICAgdmFyIHNjYWxlWSA9IDAuNTtcbiAgICAgICAgdmFyIGRlbHRhWCA9IDAuMjU7XG4gICAgICAgIHZhciBkZWx0YVkgPSAwLjI7XG5cbiAgICAgICAgLy8gc2NhbGUgcG9pbnRzIHRvIFswLDFdICsgZGVsdGFcbiAgICAgICAgdmFyIG1heFggPSBNYXRoLm1heC5hcHBseShNYXRoLCBwb2ludHMubWFwKGZ1bmN0aW9uKHApe3JldHVybiBwLnh9KSk7XG4gICAgICAgIHZhciBtaW5YID0gTWF0aC5taW4uYXBwbHkoTWF0aCwgcG9pbnRzLm1hcChmdW5jdGlvbihwKXtyZXR1cm4gcC54fSkpO1xuICAgICAgICB2YXIgbWF4WSA9IE1hdGgubWF4LmFwcGx5KE1hdGgsIHBvaW50cy5tYXAoZnVuY3Rpb24ocCl7cmV0dXJuIHAueX0pKTtcbiAgICAgICAgdmFyIG1pblkgPSBNYXRoLm1pbi5hcHBseShNYXRoLCBwb2ludHMubWFwKGZ1bmN0aW9uKHApe3JldHVybiBwLnl9KSk7XG4gICAgICAgIHBvaW50cyA9IHBvaW50cy5tYXAoZnVuY3Rpb24ocCl7XG4gICAgICAgICAgICB2YXIgeCA9IHNjYWxlWCAqIChwLngtbWluWCkvKG1heFgtbWluWCkgKyBkZWx0YVg7XG4gICAgICAgICAgICB2YXIgeSA9IHNjYWxlWSAqIChwLnktbWluWSkvKG1heFktbWluWSkgKyBkZWx0YVk7XG4gICAgICAgICAgICB2YXIgbmV3UG9pbnQgPSBuZXcgUG9pbnQoeCwgeSk7XG4gICAgICAgICAgICBuZXdQb2ludC5pZCA9IHAuaWQ7XG5cbiAgICAgICAgICAgIHJldHVybiBuZXdQb2ludDtcbiAgICAgICAgfSk7XG5cbiAgICAgICAgY29uc29sZS5sb2cocG9pbnRzKTtcblxuICAgICAgICAvLyBvbmx5IGFkZCByYW5kb20gcG9pbnRzXG4gICAgICAgIHZhciBuYlBvaW50cyA9IHBvaW50cy5sZW5ndGg7XG4gICAgICAgIGZvcih2YXIgaT0wOyBpPG5iUmFuZG9tUG9pbnRzOyArK2kpIHtcblxuICAgICAgICAgICAgdmFyIHggPSByYW5kb20oKTtcbiAgICAgICAgICAgIHZhciB5ID0gcmFuZG9tKCk7XG5cbiAgICAgICAgICAgIHZhciBuZXdQb2ludCA9IG5ldyBQb2ludCh4LCB5KTtcbiAgICAgICAgICAgIG5ld1BvaW50LmlkID0gbmJQb2ludHM7XG5cbiAgICAgICAgICAgIHBvaW50cy5wdXNoKG5ld1BvaW50KTtcblxuICAgICAgICAgICAgbmJQb2ludHMrKztcbiAgICAgICAgfVxuXG4gICAgfSBlbHNlIHtcbiAgICAgICAgLy9hZGQgcmFuZG9tIHBvaW50c1xuXG4gICAgICAgIHZhciBuYlBvaW50cyA9IDA7XG4gICAgICAgIGZvcih2YXIgaT0wOyBpPG5iUmFuZG9tUG9pbnRzOyArK2kpIHtcblxuICAgICAgICAgICAgdmFyIHggPSByYW5kb20oKTtcbiAgICAgICAgICAgIHZhciB5ID0gcmFuZG9tKCk7XG5cbiAgICAgICAgICAgIHZhciBuZXdQb2ludCA9IG5ldyBQb2ludCh4LCB5KTtcbiAgICAgICAgICAgIG5ld1BvaW50LmlkID0gbmJQb2ludHM7XG5cbiAgICAgICAgICAgIHBvaW50cy5wdXNoKG5ld1BvaW50KTtcbiAgICAgICAgICAgIFxuICAgICAgICAgICAgbmJQb2ludHMrKztcbiAgICAgICAgfVxuXG4gICAgICAgIGNpdHlTZXQgPSByYW5nZSgwLCBuYkNpdHkpO1xuICAgICAgICBjb25zb2xlLmxvZyhjaXR5U2V0KTtcbiAgICB9XG5cblxuICAgIC8vIGluaXRpYWxpemUgc3RhcnQgcG9pbnRzXG4gICAgdmFyIHBvc3NpYmxlU3RhcnRQb2ludHNJZCA9IFtdO1xuXG4gICAgZm9yICh2YXIgaSA9IDA7IGkgPCBuYlN0YXJ0UG9pbnRzOyBpKyspe1xuICAgICAgICBwb3NzaWJsZVN0YXJ0UG9pbnRzSWQucHVzaChNYXRoLmZsb29yKG5iUG9pbnRzICogcmFuZG9tKCkpKTtcbiAgICB9XG4gICAgXG5cbiAgICByZXR1cm4ge1xuICAgICAgICB0ZXh0TWVzaDogdGV4dE1lc2gsXG4gICAgICAgIHBvaW50czogcG9pbnRzLFxuICAgICAgICBjaXR5U2V0OiBjaXR5U2V0LFxuICAgICAgICBwb3NzaWJsZVN0YXJ0UG9pbnRzSWQ6IHBvc3NpYmxlU3RhcnRQb2ludHNJZCxcbiAgICAgICAgbmJSYW5kb21Qb2ludHM6IG5iUmFuZG9tUG9pbnRzLFxuICAgICAgICBmb3JjZWRFZGdlczogZm9yY2VkRWRnZXNcbiAgICB9O1xufVxuIiwiJ3VzZSBzdHJpY3QnXG5cbm1vZHVsZS5leHBvcnRzID0gZnVuY3Rpb24gKGNvbnRhaW5lcil7XG5cblx0dmFyIG1vdXNlID0ge1xuXHQgICAgeDogMCxcblx0ICAgIHk6IDBcblx0fTtcblxuXHRjb250YWluZXIuYWRkRXZlbnRMaXN0ZW5lciggJ21vdXNlbW92ZScsIGZ1bmN0aW9uKGUpe1xuXHQgICAgdmFyIHJlY3QgPSBjb250YWluZXIuZ2V0Qm91bmRpbmdDbGllbnRSZWN0KCk7XG5cblx0ICAgIG1vdXNlLnggPSAoZS5jbGllbnRYIC0gcmVjdC5sZWZ0ICkgLyByZWN0LndpZHRoO1xuXHQgICAgbW91c2UueSA9IChlLmNsaWVudFkgLSByZWN0LnRvcCApLyByZWN0LmhlaWdodDtcblx0fSk7XG5cblx0cmV0dXJuIG1vdXNlO1xuXG59O1xuIiwiJ3VzZSBzdHJpY3QnXG5cbmZ1bmN0aW9uIFBvaW50KHgsIHkpIHtcbiAgICB0aGlzLmlkID0gdW5kZWZpbmVkOyAgICAgICAgICAgICAgICBcbiAgICB0aGlzLnggPSB4O1xuICAgIHRoaXMueSA9IHk7XG4gICAgdGhpcy5uZXh0cyA9IFtdO1xufVxuXG5tb2R1bGUuZXhwb3J0cyA9IFBvaW50OyIsIid1c2Ugc3RyaWN0J1xuXG52YXIgYW50RnVuY3Rpb24gPSByZXF1aXJlKCcuL2FudC5qcycpO1xudmFyIGFudHNHcm91cCA9IHJlcXVpcmUoJy4vYW50c0dyb3VwJyk7XG5cbnZhciByYW5kb20gPSBNYXRoLnJhbmRvbTtcblxudmFyIFJBTkRPTU1WVCA9IDAuMDAzO1xudmFyIEFOVFNJWkUgPSAwLjAwMjtcblxubW9kdWxlLmV4cG9ydHMgPSBmdW5jdGlvbihjb250YWluZXIsIHBvaW50c01hcCwgb3B0aW9ucyl7XG5cblx0aWYoIWNvbnRhaW5lcilcblx0XHR0aHJvdyBuZXcgVHlwZUVycm9yKCdNaXNzaW5nIGNvbnRhaW5lcicpO1xuXG5cdC8vIEFudHMgdmFyaWFibGVzXG5cdHZhciBlZGdlcyA9IHBvaW50c01hcC5lZGdlcztcblx0dmFyIG9ialBvcHVsYXRpb25Jbml0aWFsID0gb3B0aW9ucy5uYkFudHM7XG5cdHZhciBvYmpQb3B1bGF0aW9uID0gb2JqUG9wdWxhdGlvbkluaXRpYWw7XG5cdHZhciBwb2ludHNJbmZvcyA9IHBvaW50c01hcC5wb2ludHNJbmZvcztcblx0dmFyIHBvcHVsYXRpb24gPSBbXTtcblx0dmFyIG5iQW50c1BlclN0ZXAgPSAxMDA7XG5cdFxuXHR2YXIgQW50ID0gYW50RnVuY3Rpb24oY29udGFpbmVyLCBwb2ludHNJbmZvcywgb3B0aW9ucyk7XG5cdGFudHNHcm91cCA9IGFudHNHcm91cChBbnQpO1xuXG5cdC8vIEFuaW1hdGlvbiB2YXJpYWJsZXNcblx0dmFyIGFuaW1JRDtcblx0dmFyIGRlbHRhVGltZTtcblx0dmFyIEZQU0NvdW50O1xuXHR2YXIgbGFzdFVwZGF0ZSA9IHBlcmZvcm1hbmNlLm5vdygpO1xuXHR2YXIgRlBTTW9uaXRvciA9IGRvY3VtZW50LnF1ZXJ5U2VsZWN0b3IoJyNGUFMnKTtcblx0dmFyIGRUTW9uaXRvciA9IGRvY3VtZW50LnF1ZXJ5U2VsZWN0b3IoJyNkVCcpO1xuXHR2YXIgcmVmcmVzaFRpbWUgPSAwO1xuXHR2YXIgbWF4RGVsdGFUaW1lID0gMzA7XG5cdHZhciBGUFNPdmVyTGltaXRDb3VudCA9IDA7XG5cdHZhciBGUFNVbmRlckxpbWl0Q291bnQgPSAwO1xuXG5cblx0Ly8gQ2FudmFzXG5cdHZhciBjYW52YXNMaXN0ID0gZG9jdW1lbnQuZ2V0RWxlbWVudHNCeVRhZ05hbWUoXCJjYW52YXNcIik7XG5cdFxuXHRpZiAoY2FudmFzTGlzdC5sZW5ndGggPT09IDApe1xuXHRcdHZhciBjYW52YXMgPSBkb2N1bWVudC5jcmVhdGVFbGVtZW50KFwiY2FudmFzXCIpO1xuXHRcdHZhciByZWN0ID0gY29udGFpbmVyLmdldEJvdW5kaW5nQ2xpZW50UmVjdCgpO1xuXHRcdGNhbnZhcy53aWR0aCA9IHJlY3Qud2lkdGg7XG5cdFx0Y2FudmFzLmhlaWdodCA9IHJlY3QuaGVpZ2h0O1xuXHRcdGNhbnZhcy5zdHlsZS5iYWNrZ3JvdW5kQ29sb3IgPSBcInJnYmEoMjUwLCAyNTAsIDI1MCwgMClcIjsgXG5cdFx0Y29udGFpbmVyLmFwcGVuZENoaWxkKGNhbnZhcyk7XG5cdH1cblx0ZWxzZXtcblx0XHR2YXIgY2FudmFzID0gY2FudmFzTGlzdFswXTtcblx0XHRjb25zb2xlLmxvZygnQ0FOVkFTJyk7XG5cdH1cblx0dmFyIGNvbnRleHQgPSBjYW52YXMuZ2V0Q29udGV4dChcIjJkXCIpO1xuXHRjb250ZXh0LmNsZWFyUmVjdCAoIDAgLCAwICwgY2FudmFzLndpZHRoLCBjYW52YXMuaGVpZ2h0ICk7XG5cdFxuXG5cdGZ1bmN0aW9uIGNoZWNrQW50TnVtYmVyKGFudE51bWJlcil7XG5cdFx0aWYgKGFudE51bWJlciA8IG9ialBvcHVsYXRpb24gLSA1MCl7XG5cdFx0XHRGUFNNb25pdG9yLnN0eWxlLmNvbG9yID0gXCJncmVlblwiO1xuXHRcdFx0cG9wdWxhdGlvbiA9IGFudHNHcm91cC5jcmVhdGUocG9wdWxhdGlvbik7XG5cdFx0fVx0XG5cdFx0ZWxzZSBpZiAoYW50TnVtYmVyID4gb2JqUG9wdWxhdGlvbil7XG5cdFx0XHRwb3B1bGF0aW9uID0gYW50c0dyb3VwLnJlbW92ZShwb3B1bGF0aW9uLCBhbnROdW1iZXIgLSBvYmpQb3B1bGF0aW9uKTtcblx0XHRcdEZQU01vbml0b3Iuc3R5bGUuY29sb3IgPSBcInJlZFwiO1xuXHRcdH1cblx0XHRlbHNlXG5cdFx0XHRGUFNNb25pdG9yLnN0eWxlLmNvbG9yID0gXCJ3aGl0ZVwiO1xuXHR9XG5cblx0ZnVuY3Rpb24gZGlzcGxheUZQUyhkVCl7XG5cdFx0RlBTQ291bnQgPSAoMTAwMC9kVCkudG9GaXhlZCgyKTtcblx0XHR2YXIgdCA9IGRULnRvRml4ZWQoMik7XG5cdFx0RlBTTW9uaXRvci5pbm5lclRleHQgPSAnRlBTIDogJyArIEZQU0NvdW50OyAgXG5cdFx0ZFRNb25pdG9yLmlubmVyVGV4dCA9ICduYkFudHMgOiAnICsgcG9wdWxhdGlvbi5sZW5ndGg7XG5cdFx0Ly8gZFRNb25pdG9yLmlubmVyVGV4dCA9ICdkVCA6ICcgKyB0ICsgJ21zJztcblx0fVxuXG5cdGZ1bmN0aW9uIHRpY2soKSB7XG5cdFx0dmFyIG5vdyA9IHBlcmZvcm1hbmNlLm5vdygpO1xuXHRcdGRlbHRhVGltZSA9IG5vdyAtIGxhc3RVcGRhdGU7XG5cdFx0bGFzdFVwZGF0ZSA9IG5vdztcblx0XHRyZWZyZXNoVGltZSArPSBkZWx0YVRpbWUvMTAwMDsgLy8gaW4gc2Vjb25kc1xuXG5cdFx0Y29uc29sZS5sb2coJ25iQW50cycsIHBvcHVsYXRpb24ubGVuZ3RoKTtcblxuXHRcdGNoZWNrQW50TnVtYmVyKHBvcHVsYXRpb24ubGVuZ3RoKTtcblxuXHRcdC8vIGRpc3BsYXkgRlBTIGluZm8gZXZlcnkgMC4zIHNcblx0XHRpZiAocmVmcmVzaFRpbWUgPiAwLjMpe1xuXHRcdFx0ZGlzcGxheUZQUyhkZWx0YVRpbWUpO1xuXHRcdFx0cmVmcmVzaFRpbWUgPSAwOyBcblx0XHR9XG5cblx0XHQvLyByZW1vdmUgYW50cyB3aGVuIGZyYW1lIHJhdGUgaXMgdG9vIGxvd1xuXHRcdGlmIChGUFNPdmVyTGltaXRDb3VudCA9PT0gMTApIHtcblx0XHRcdG9ialBvcHVsYXRpb24gPSBvYmpQb3B1bGF0aW9uICogbWF4RGVsdGFUaW1lIC8gZGVsdGFUaW1lO1xuXHRcdFx0RlBTT3ZlckxpbWl0Q291bnQgPSAwO1xuXHRcdH1cblxuXHRcdHdoaWxlIChGUFNVbmRlckxpbWl0Q291bnQgPiA1MCAmJiBvYmpQb3B1bGF0aW9uIDwgb2JqUG9wdWxhdGlvbkluaXRpYWwpIHtcblx0XHRcdG9ialBvcHVsYXRpb24gKz0gMTAwO1xuXHRcdFx0Ly8gRlBTT3ZlckxpbWl0Q291bnQgPSAwO1xuXHRcdH1cblxuXHRcdC8vIGNoZWNrIGR1cmF0aW9uIG9mIG92ZXIvdW5kZXIgZnJhbWVyYXRlIGxpbWl0IHBlcmlvZHNcblx0XHRpZiAoZGVsdGFUaW1lID4gbWF4RGVsdGFUaW1lKXtcblx0XHRcdEZQU092ZXJMaW1pdENvdW50Kys7XG5cdFx0XHRGUFNVbmRlckxpbWl0Q291bnQgPSAwO1xuXHRcdH1cblx0XHRlbHNlIHtcblx0XHRcdEZQU092ZXJMaW1pdENvdW50ID0gMDtcblx0XHRcdEZQU1VuZGVyTGltaXRDb3VudCsrO1xuXHRcdH1cblxuXHRcdC8vIGRyYXcgaW4gY2FudmFzXG5cdFx0dmFyIHcgPSBjYW52YXMud2lkdGg7XG5cdFx0dmFyIGggPSBjYW52YXMuaGVpZ2h0O1xuXHRcdHZhciBtb3VzZSA9IFtsYXN0TW91c2VNb3ZlRXZlbnQuY2xpZW50WC93LCBsYXN0TW91c2VNb3ZlRXZlbnQuY2xpZW50WS9oXTtcblx0XHRjb250ZXh0LnNldFRyYW5zZm9ybSh3LCAwLCAwLCBoLCAwLCAwKTtcblx0XHRjb250ZXh0LmZpbGxTdHlsZSA9IFwicmdiYSgyNTAsIDI1MCwgMjUwLCAwLjQpXCI7XG5cdFx0Y29udGV4dC5maWxsUmVjdCgwLDAsdyxoKTtcblxuXHRcdC8vIGVkZ2VzXG5cdFx0Ly8gY29udGV4dC5zdHJva2VTdHlsZSA9IFwiIzAwMFwiO1xuXHRcdC8vIGZvcih2YXIgaT0wOyBpIDwgZWRnZXMubGVuZ3RoOyArK2kpIHtcblx0XHQvLyAgICAgY29udGV4dC5saW5lV2lkdGggPSAwLjAwMDE7XG5cdFx0Ly8gICAgIHZhciBlZGdlID0gZWRnZXNbaV07XG5cdFx0Ly8gICAgIC8vIGlmIChlZGdlLnBoZXJvbW9uICE9IDApe1xuXHRcdC8vICAgICAvLyAgICAgY29udGV4dC5saW5lV2lkdGggPSBNYXRoLm1pbigwLjAwMDAxICogZWRnZS5waGVyb21vbiwgMC4wMSk7XG5cdFx0Ly8gICAgIC8vIH0gZWxzZSB7XG5cdFx0Ly8gICAgIC8vICAgICBjb250ZXh0LmxpbmVXaWR0aCA9IDAuMDAwMDE7XG5cdFx0Ly8gICAgIC8vIH1cblx0XHQvLyAgICAgY29udGV4dC5iZWdpblBhdGgoKTtcblx0XHQvLyAgICAgY29udGV4dC5tb3ZlVG8ocG9pbnRzW2VkZ2UucHQxLmlkXS54LCBwb2ludHNbZWRnZS5wdDEuaWRdLnkpO1xuXHRcdC8vICAgICBjb250ZXh0LmxpbmVUbyhwb2ludHNbZWRnZS5wdDIuaWRdLngsIHBvaW50c1tlZGdlLnB0Mi5pZF0ueSk7XG5cdFx0Ly8gICAgIGNvbnRleHQuc3Ryb2tlKCk7XG5cdFx0Ly8gfVxuXG5cdFx0Ly8gLy8gdmVydGljZXNcblx0XHQvLyBmb3IodmFyIGk9MDsgaTxwb2ludHMubGVuZ3RoOyArK2kpIHtcblx0XHQvLyAgICAgY29udGV4dC5iZWdpblBhdGgoKVxuXHRcdC8vICAgICB2YXIgcG9pbnQgPSBwb2ludHNbaV07XG5cdFx0Ly8gICAgIGlmIChjaXR5U2V0Lmhhcyhwb2ludC5pZCkpIHtcblx0XHQvLyAgICAgICAgIGNvbnRleHQuZmlsbFN0eWxlID0gXCIjMDEwMURGXCI7XG5cdFx0Ly8gICAgICAgICBjb250ZXh0LmFyYyhwb2ludC54LCBwb2ludC55LCAwLjAwNiwgMCwgMipNYXRoLlBJKTtcblx0XHQvLyAgICAgfVxuXHRcdC8vICAgICBlbHNlIHtcblx0XHQvLyAgICAgICAgIGNvbnRleHQuZmlsbFN0eWxlID0gXCIjMDAwXCI7XG5cdFx0Ly8gICAgICAgICBjb250ZXh0LmFyYyhwb2ludHNbaV0ueCwgcG9pbnRzW2ldLnksIDAuMDAzLCAwLCAyKk1hdGguUEkpO1xuXHRcdC8vICAgICB9XG5cdFx0Ly8gICAgIGNvbnRleHQuY2xvc2VQYXRoKCk7XG5cdFx0Ly8gICAgIGNvbnRleHQuZmlsbCgpO1xuXHRcdC8vIH1cblxuXHRcdC8vIG1vdmUgYW50c1xuXHRcdHBvcHVsYXRpb24uZm9yRWFjaChmdW5jdGlvbihhbnQpe1xuXHRcdFx0YW50LnRyYW5zaXQoKTtcblx0XHR9KTtcblxuXHRcdC8vIHBoZXJvbW9uIGV2YXBvcmF0aW9uXG5cdFx0ZWRnZXMuZm9yRWFjaChmdW5jdGlvbihlZGdlKXtcblx0XHRcdGlmKGVkZ2UucGhlcm9tb24gPiAwKXtcblx0XHRcdFx0ZWRnZS5waGVyb21vbiAtPSAwLjAwMDE7XG5cdFx0XHR9XG5cdFx0fSk7XG5cblx0XHQvLyBhbnRzXG5cdFx0cG9wdWxhdGlvbi5mb3JFYWNoKGZ1bmN0aW9uKGFudCl7XG5cdFx0XHRjb250ZXh0LmJlZ2luUGF0aCgpXG5cdFx0XHR2YXIgeCA9IGFudC54ICsgUkFORE9NTVZUKnJhbmRvbSgpO1xuXHRcdFx0dmFyIHkgPSBhbnQueSArIFJBTkRPTU1WVCpyYW5kb20oKTtcblxuXHRcdFx0Y29udGV4dC5maWxsU3R5bGUgPSBcImJsYWNrXCJcblx0XHRcdGNvbnRleHQuZmlsbFJlY3QoeCwgeSwgQU5UU0laRSwgQU5UU0laRSk7XG5cdFx0XHRjb250ZXh0LmNsb3NlUGF0aCgpO1xuXHRcdFx0Y29udGV4dC5maWxsKCk7XG5cdFx0fSlcblx0fTtcblx0XG5cdHZhciBsYXN0TW91c2VNb3ZlRXZlbnQgPSB7XG5cdFx0Y2xpZW50WDogMCxcblx0XHRjbGllbnRZOiAwXG5cdH07XG5cdFxuXHRjb250YWluZXIuYWRkRXZlbnRMaXN0ZW5lcignbW91c2Vtb3ZlJywgZnVuY3Rpb24oZSl7XG5cdFx0bGFzdE1vdXNlTW92ZUV2ZW50ID0gZTtcblx0fSk7XG5cdFxuXHR2YXIgcGF1c2VkID0gZmFsc2U7XG5cdFxuXHRmdW5jdGlvbiB0b2dnbGVQbGF5UGF1c2UoKXtcblx0XHRwYXVzZWQgPSAhcGF1c2VkO1xuXHRcdGlmKCFwYXVzZWQpXG5cdFx0XHRhbmltYXRlKCk7XG5cdH1cblxuXHRmdW5jdGlvbiByZXNldCgpe1xuXHRcdHBvcHVsYXRpb24gPSBbXTtcblx0XHRlZGdlcyA9IFtdO1xuXHRcdHBvaW50c0luZm9zID0gW107XG5cblx0XHRjYW5jZWxBbmltYXRpb25GcmFtZShhbmltSUQpO1xuXHR9XG5cdFxuXHQvLyBjb250YWluZXIuYWRkRXZlbnRMaXN0ZW5lcignY2xpY2snLCB0b2dnbGVQbGF5UGF1c2UpO1xuXG5cdGZ1bmN0aW9uIGFuaW1hdGUoKXtcblx0XHR0aWNrKCk7XG5cdFx0XG5cdFx0aWYoIXBhdXNlZClcblx0XHRcdGFuaW1JRCA9IHJlcXVlc3RBbmltYXRpb25GcmFtZShhbmltYXRlKTtcblx0fVxuXHRhbmltYXRlKCk7XG5cblx0ZnVuY3Rpb24gbW9kaWZ5QW50cyhvcHRzKXtcblx0XHRvYmpQb3B1bGF0aW9uID0gb3B0cy5uYkFudHM7XG5cblx0XHRwb3B1bGF0aW9uLmZvckVhY2goZnVuY3Rpb24oYW50KXtcblx0XHRcdGFudC52ZWxvY2l0eSA9IG9wdHMudmVsb2NpdHk7XG5cdFx0XHRhbnQud2VpZ2h0ID0gb3B0cy53ZWlnaHQ7XG5cdFx0XHRhbnQucmVwU2l6ZSA9IG9wdHMucmVwU2l6ZTtcblx0XHRcdGFudC5yZXBTcGVlZCA9IG9wdHMucmVwU3BlZWQ7XG5cdFx0fSk7XG5cdH1cblx0XG5cdHJldHVybiB7XG5cdFx0dG9nZ2xlUGxheVBhdXNlOiB0b2dnbGVQbGF5UGF1c2UsXG5cdFx0cmVzZXQ6IHJlc2V0LFxuXHRcdC8vIHNob3VsZCBiZSBhIGdldHRlci9zZXR0ZXIsIGJ1dCBJRThcblx0XHRnZXRBbnRDb3VudDogZnVuY3Rpb24oKXtcblx0XHRcdHJldHVybiBwb3B1bGF0aW9uLmxlbmd0aDtcblx0XHR9LFxuXHRcdG1vZGlmeUFudHM6IG1vZGlmeUFudHNcblx0fVxufVxuIiwiJ3VzZSBzdHJpY3QnO1xuXG52YXIgcG9pbnQxID0gXCJtIDE4LjI1LDE5LjUgYyAxOC4yNSwxNS43Njc3Njg4IDE1LjAwMDA5NTEsMTIuNzUgMTEsMTIuNzUgYyA2Ljk5OTkwNDg4LDEyLjc1IDMuNzUsMTUuNzY3NzY4OCAzLjc1LDE5LjUgYyAzLjc1LDIzLjIzMjIzMTIgNi45OTk5MDQ4OCwyNi4yNSAxMSwyNi4yNSBjIDE1LjAwMDA5NTEsMjYuMjUgMTguMjUsMjMuMjMyMjMxMiAxOC4yNSwxOS41IHogbSA0LjI1LDE5LjUgYyA0LjI1LDE2LjA1MjUyOTQgNy4yNjgxMDg2MywxMy4yNSAxMSwxMy4yNSBjIDE0LjczMTg5MTQsMTMuMjUgMTcuNzUsMTYuMDUyNTI5NCAxNy43NSwxOS41IGMgMTcuNzUsMjIuOTQ3NDcwNiAxNC43MzE4OTE0LDI1Ljc1IDExLDI1Ljc1IGMgNy4yNjgxMDg2MywyNS43NSA0LjI1LDIyLjk0NzQ3MDYgNC4yNSwxOS41IHpcIjtcbnZhciBwb2ludDIgPSBcIm0gODkuMjUsOC41IGMgODkuMjUsNC4yMTYwNTc4NyA4NS41NTI4NzE2LDAuNzUgODEsMC43NSBjIDc2LjQ0NzEyODQsMC43NSA3Mi43NSw0LjIxNjA1Nzg3IDcyLjc1LDguNSBjIDcyLjc1LDEyLjc4Mzk0MjEgNzYuNDQ3MTI4NCwxNi4yNSA4MSwxNi4yNSBjIDg1LjU1Mjg3MTYsMTYuMjUgODkuMjUsMTIuNzgzOTQyMSA4OS4yNSw4LjUgeiBtIDczLjI1LDguNSBjIDczLjI1LDQuNDk5NjcwODggNzYuNzE2MzE1NiwxLjI1IDgxLDEuMjUgYyA4NS4yODM2ODQ0LDEuMjUgODguNzUsNC40OTk2NzA4OCA4OC43NSw4LjUgYyA4OC43NSwxMi41MDAzMjkxIDg1LjI4MzY4NDQsMTUuNzUgODEsMTUuNzUgYyA3Ni43MTYzMTU2LDE1Ljc1IDczLjI1LDEyLjUwMDMyOTEgNzMuMjUsOC41IHpcIjtcbnZhciBwb2ludDMgPSBcIm0gMTYwLjI1LDExIGMgMTYwLjI1LDUuMzMzNDgzNzUgMTU1LjIwODE2OCwwLjc1IDE0OSwwLjc1IGMgMTQyLjc5MTgzMiwwLjc1IDEzNy43NSw1LjMzMzQ4Mzc1IDEzNy43NSwxMSBjIDEzNy43NSwxNi42NjY1MTYyIDE0Mi43OTE4MzIsMjEuMjUgMTQ5LDIxLjI1IGMgMTU1LjIwODE2OCwyMS4yNSAxNjAuMjUsMTYuNjY2NTE2MiAxNjAuMjUsMTEgeiBtIDEzOC4yNSwxMSBjIDEzOC4yNSw1LjYyMDgyMTI1IDE0My4wNTc5MDMsMS4yNSAxNDksMS4yNSBjIDE1NC45NDIwOTcsMS4yNSAxNTkuNzUsNS42MjA4MjEyNSAxNTkuNzUsMTEgYyAxNTkuNzUsMTYuMzc5MTc4NyAxNTQuOTQyMDk3LDIwLjc1IDE0OSwyMC43NSBjIDE0My4wNTc5MDMsMjAuNzUgMTM4LjI1LDE2LjM3OTE3ODcgMTM4LjI1LDExIHpcIjtcbnZhciBwb2ludDQgPSBcIm0gMTYwLjI1LDc2LjUgYyAxNjAuMjUsNzIuNzY3NzY4OCAxNTcuMDAwMDk1LDY5Ljc1IDE1Myw2OS43NSBjIDE0OC45OTk5MDUsNjkuNzUgMTQ1Ljc1LDcyLjc2Nzc2ODggMTQ1Ljc1LDc2LjUgYyAxNDUuNzUsODAuMjMyMjMxMiAxNDguOTk5OTA1LDgzLjI1IDE1Myw4My4yNSBjIDE1Ny4wMDAwOTUsODMuMjUgMTYwLjI1LDgwLjIzMjIzMTIgMTYwLjI1LDc2LjUgeiBtIDE0Ni4yNSw3Ni41IGMgMTQ2LjI1LDczLjA1MjUyOTQgMTQ5LjI2ODEwOSw3MC4yNSAxNTMsNzAuMjUgYyAxNTYuNzMxODkxLDcwLjI1IDE1OS43NSw3My4wNTI1Mjk0IDE1OS43NSw3Ni41IGMgMTU5Ljc1LDc5Ljk0NzQ3MDYgMTU2LjczMTg5MSw4Mi43NSAxNTMsODIuNzUgYyAxNDkuMjY4MTA5LDgyLjc1IDE0Ni4yNSw3OS45NDc0NzA2IDE0Ni4yNSw3Ni41IHpcIjtcbnZhciBwb2ludDUgPSBcIm0gOTUuMjUsNzYgYyA5NS4yNSw3MC4zMzYyNzk1IDkwLjQzNDQwNjUsNjUuNzUgODQuNSw2NS43NSBjIDc4LjU2NTU5MzUsNjUuNzUgNzMuNzUsNzAuMzM2Mjc5NSA3My43NSw3NiBjIDczLjc1LDgxLjY2MzcyMDUgNzguNTY1NTkzNSw4Ni4yNSA4NC41LDg2LjI1IGMgOTAuNDM0NDA2NSw4Ni4yNSA5NS4yNSw4MS42NjM3MjA1IDk1LjI1LDc2IHogbSA3NC4yNSw3NiBjIDc0LjI1LDcwLjYxODAyNTUgNzguODM2NDI2OCw2Ni4yNSA4NC41LDY2LjI1IGMgOTAuMTYzNTczMiw2Ni4yNSA5NC43NSw3MC42MTgwMjU1IDk0Ljc1LDc2IGMgOTQuNzUsODEuMzgxOTc0NSA5MC4xNjM1NzMyLDg1Ljc1IDg0LjUsODUuNzUgYyA3OC44MzY0MjY4LDg1Ljc1IDc0LjI1LDgxLjM4MTk3NDUgNzQuMjUsNzYgelwiO1xudmFyIHBvaW50NiA9IFwibSAyMC4yNSw3NSBjIDIwLjI1LDcwLjk5MTkzMzggMTYuNzc2NDk5NSw2Ny43NSAxMi41LDY3Ljc1IGMgOC4yMjM1MDA0Niw2Ny43NSA0Ljc1LDcwLjk5MTkzMzggNC43NSw3NSBjIDQuNzUsNzkuMDA4MDY2MiA4LjIyMzUwMDQ2LDgyLjI1IDEyLjUsODIuMjUgYyAxNi43NzY0OTk1LDgyLjI1IDIwLjI1LDc5LjAwODA2NjIgMjAuMjUsNzUgeiBtIDUuMjUsNzUgYyA1LjI1LDcxLjI3NjA3OTcgOC40OTIyMjgyOSw2OC4yNSAxMi41LDY4LjI1IGMgMTYuNTA3NzcxNyw2OC4yNSAxOS43NSw3MS4yNzYwNzk3IDE5Ljc1LDc1IGMgMTkuNzUsNzguNzIzOTIwMyAxNi41MDc3NzE3LDgxLjc1IDEyLjUsODEuNzUgYyA4LjQ5MjIyODI5LDgxLjc1IDUuMjUsNzguNzIzOTIwMyA1LjI1LDc1IHpcIjtcbnZhciBwb2ludDcgPSBcIm0gMjMuMjUsMTM5IGMgMjMuMjUsMTMyLjc4Njc5NyAxOC4yMTMyMDM0LDEyNy43NSAxMiwxMjcuNzUgYyA1Ljc4Njc5NjU2LDEyNy43NSAwLjc1LDEzMi43ODY3OTcgMC43NSwxMzkgYyAwLjc1LDE0NS4yMTMyMDMgNS43ODY3OTY1NiwxNTAuMjUgMTIsMTUwLjI1IGMgMTguMjEzMjAzNCwxNTAuMjUgMjMuMjUsMTQ1LjIxMzIwMyAyMy4yNSwxMzkgeiBtIDEuMjUsMTM5IGMgMS4yNSwxMzMuMDYyOTM5IDYuMDYyOTM4OTQsMTI4LjI1IDEyLDEyOC4yNSBjIDE3LjkzNzA2MTEsMTI4LjI1IDIyLjc1LDEzMy4wNjI5MzkgMjIuNzUsMTM5IGMgMjIuNzUsMTQ0LjkzNzA2MSAxNy45MzcwNjExLDE0OS43NSAxMiwxNDkuNzUgYyA2LjA2MjkzODk0LDE0OS43NSAxLjI1LDE0NC45MzcwNjEgMS4yNSwxMzkgelwiO1xudmFyIHBvaW50OCA9IFwibSA5NS4yNDI5MjA5LDEzMy4xNDE2MjQgYyA5NS41NTM3OTI3LDEyNy4yMDk4MzYgOTAuNzU2MjY5OCwxMjIuMTQxNzM2IDg0LjUzMzc0MDQsMTIxLjgxNTYyNyBjIDc4LjMxMTIxMSwxMjEuNDg5NTE4IDczLjAxMDIwODYsMTI2LjAyODM3NyA3Mi42OTkzMzY4LDEzMS45NjAxNjUgYyA3Mi4zODg0NjUsMTM3Ljg5MTk1MyA3Ny4xODU5ODc5LDE0Mi45NjAwNTMgODMuNDA4NTE3MywxNDMuMjg2MTYyIGMgODkuNjMxMDQ2NywxNDMuNjEyMjcxIDk0LjkzMjA0OSwxMzkuMDczNDEyIDk1LjI0MjkyMDksMTMzLjE0MTYyNCB6IG0gNzMuMTk4NjUxNiwxMzEuOTg2MzMzIGMgNzMuNDk0NzcxMSwxMjYuMzM2MDM2IDc4LjU1NTM4OCwxMjIuMDAzMDAxIDg0LjUwNzU3MjQsMTIyLjMxNDk0MiBjIDkwLjQ1OTc1NjcsMTIyLjYyNjg4MyA5NS4wMzk3MjU2LDEyNy40NjUxNTkgOTQuNzQzNjA2MSwxMzMuMTE1NDU2IGMgOTQuNDQ3NDg2NiwxMzguNzY1NzU0IDg5LjM4Njg2OTYsMTQzLjA5ODc4OCA4My40MzQ2ODUzLDE0Mi43ODY4NDcgYyA3Ny40ODI1MDA5LDE0Mi40NzQ5MDcgNzIuOTAyNTMyLDEzNy42MzY2MyA3My4xOTg2NTE2LDEzMS45ODYzMzMgelwiO1xudmFyIHBvaW50OSA9IFwibSAxNjcuNzI4MDIzLDEzNS4wNjI0NzkgYyAxNjguMDM4MzI3LDEyOS4xNDE1MzMgMTYzLjQ4MjExNywxMjQuMDg5ODY0IDE1Ny41NTE2NjIsMTIzLjc3OTA2MiBjIDE1MS42MjEyMDgsMTIzLjQ2ODI2MSAxNDYuNTYxOTE0LDEyOC4wMTYwMDIgMTQ2LjI1MTYxLDEzMy45MzY5NDggYyAxNDUuOTQxMzA3LDEzOS44NTc4OTQgMTUwLjQ5NzUxNywxNDQuOTA5NTYyIDE1Ni40Mjc5NzEsMTQ1LjIyMDM2NCBjIDE2Mi4zNTg0MjYsMTQ1LjUzMTE2NiAxNjcuNDE3NzIsMTQwLjk4MzQyNSAxNjcuNzI4MDIzLDEzNS4wNjI0NzkgeiBtIDE0Ni43NTA5MjUsMTMzLjk2MzExNiBjIDE0Ny4wNDY3NjcsMTI4LjMxODEyMSAxNTEuODcwNjE3LDEyMy45ODIwMTggMTU3LjUyNTQ5NCwxMjQuMjc4Mzc3IGMgMTYzLjE4MDM3MiwxMjQuNTc0NzM3IDE2Ny41MjQ1NSwxMjkuMzkxMzE3IDE2Ny4yMjg3MDksMTM1LjAzNjMxMSBjIDE2Ni45MzI4NjcsMTQwLjY4MTMwNSAxNjIuMTA5MDE3LDE0NS4wMTc0MDkgMTU2LjQ1NDEzOSwxNDQuNzIxMDUgYyAxNTAuNzk5MjYyLDE0NC40MjQ2OSAxNDYuNDU1MDgzLDEzOS42MDgxMSAxNDYuNzUwOTI1LDEzMy45NjMxMTYgelwiO1xudmFyIGxldHRyZXNBTlQgPSBcIm0gMjUuNjE4NTQxNCwzMS45NDQ3MjY2IGwgMjUuNjE4NTQxNCwzMS4xMDAzOTA2IGwgMjUuMzg4MjY4LDMwLjQwOTU3MDMgbCAyNC45Mjc3MjExLDMwLjAyNTc4MTMgbCAyNC41NDM5MzIxLDI5LjU2NTIzNDQgbCAyNC4wODMzODUyLDI5LjE4MTQ0NTMgbCAyMy42OTk1OTYxLDI4LjcyMDg5ODQgbCAyMy4wMDg3NzU4LDI4LjQ5MDYyNSBsIDIwLjQ3NTc2OCwyOC40OTA2MjUgbCAxOS44NjE3MDU1LDI4LjcyMDg5ODQgbCAxOS40Nzc5MTY0LDI5LjE4MTQ0NTMgbCAxOC43ODcwOTYxLDI5LjMzNDk2MDkgbCAxNC41NjU0MTY0LDI5LjMzNDk2MDkgbCAxMy45NTEzNTM5LDI5LjU2NTIzNDQgbCAxMy41Njc1NjQ5LDMwLjAyNTc4MTMgbCAxMi44NzY3NDQ2LDMwLjE3OTI5NjkgbCAxMi4wMzI0MDg2LDMwLjE3OTI5NjkgbCAxMS40MTgzNDYxLDMwLjQwOTU3MDMgbCAxMS4wMzQ1NTcxLDMwLjg3MDExNzIgbCAxMC4zNDM3MzY3LDMxLjAyMzYzMjggbCA5LjcyOTY3NDI0LDMxLjI1MzkwNjMgbCA5LjM0NTg4NTE4LDMxLjcxNDQ1MzEgbCA4Ljg4NTMzODMsMzIuMDk4MjQyMiBsIDguNTAxNTQ5MjQsMzIuNTU4Nzg5MSBsIDguMDQxMDAyMzcsMzIuOTQyNTc4MSBsIDcuODg3NDg2NzQsMzMuNjMzMzk4NCBsIDcuNjU3MjEzMywzNC4yNDc0NjA5IGwgNy4xOTY2NjY0MywzNC42MzEyNSBsIDcuMDQzMTUwOCwzNS4zMjIwNzAzIGwgNy4wNDMxNTA4LDM2LjE2NjQwNjMgbCA3LjE5NjY2NjQzLDM2Ljc4MDQ2ODggbCA3LjY1NzIxMzMsMzcuMTY0MjU3OCBsIDguMDQxMDAyMzcsMzcuNjI0ODA0NyBsIDguNTAxNTQ5MjQsMzguMDA4NTkzOCBsIDguODg1MzM4MywzOC40NjkxNDA2IGwgOS4zNDU4ODUxOCwzOC44NTI5Mjk3IGwgOS43Mjk2NzQyNCwzOS4zMTM0NzY2IGwgMTAuMzQzNzM2NywzOS40NjY5OTIyIGwgMTUuNDA5NzUyNCwzOS40NjY5OTIyIGwgMTYuMTAwNTcyNywzOS4zMTM0NzY2IGwgMTYuNDg0MzYxNywzOC44NTI5Mjk3IGwgMTcuMDk4NDI0MiwzOC42MjI2NTYzIGwgMTcuNzg5MjQ0NiwzOC40NjkxNDA2IGwgMTguMTczMDMzNiwzOC4wMDg1OTM4IGwgMTguNjMzNTgwNSwzNy42MjQ4MDQ3IGwgMTkuMDE3MzY5NiwzNy4xNjQyNTc4IGwgMTkuNjMxNDMyMSwzNi45MzM5ODQ0IGwgMjAuMzIyMjUyNCwzNi43ODA0Njg4IGwgMjAuNzA2MDQxNCwzNi4zMTk5MjE5IGwgMjEuMzIwMTAzOSwzNi4wODk2NDg0IGwgMjIuMDEwOTI0MiwzNS45MzYxMzI4IGwgMjIuMzk0NzEzMywzNS40NzU1ODU5IGwgMjMuMDA4Nzc1OCwzNS4yNDUzMTI1IGwgMjMuNjk5NTk2MSwzNS4wOTE3OTY5IGwgMjMuOTI5ODY5NiwzNC40Nzc3MzQ0IGwgMjQuMDgzMzg1MiwzMy43ODY5MTQxIGwgMjQuNTQzOTMyMSwzMy40MDMxMjUgbCAyNC45Mjc3MjExLDMyLjk0MjU3ODEgbCAyNS4zODgyNjgsMzIuNTU4Nzg5MSBsIDI1LjYxODU0MTQsMzEuOTQ0NzI2NiAgbSAzNy41OTI3NjAyLDQzLjkxODk0NTMgbCAzNy4yMDg5NzExLDQ0LjM3OTQ5MjIgbCAzNi43NDg0MjQyLDQ0Ljc2MzI4MTMgbCAzNi4zNjQ2MzUyLDQ1LjIyMzgyODEgbCAzNS42NzM4MTQ5LDQ1LjM3NzM0MzggbCAzNS4wNTk3NTI0LDQ1LjYwNzYxNzIgbCAzNC42NzU5NjMzLDQ2LjA2ODE2NDEgbCAzMy45ODUxNDMsNDYuMjIxNjc5NyBsIDMxLjQ1MjEzNTIsNDYuMjIxNjc5NyBsIDMwLjgzODA3MjcsNDYuMDY4MTY0MSBsIDMwLjQ1NDI4MzYsNDUuNjA3NjE3MiBsIDI5Ljc2MzQ2MzMsNDUuMzc3MzQzOCBsIDI5LjE0OTQwMDgsNDUuMjIzODI4MSBsIDI4Ljc2NTYxMTcsNDQuNzYzMjgxMyBsIDI4LjMwNTA2NDksNDQuMzc5NDkyMiBsIDI3LjkyMTI3NTgsNDMuOTE4OTQ1MyBsIDI3LjIzMDQ1NTUsNDMuNjg4NjcxOSBsIDI2LjYxNjM5Myw0My41MzUxNTYzIGwgMjYuMjMyNjAzOSw0My4wNzQ2MDk0IGwgMjUuNzcyMDU3MSw0My4wNzQ2MDk0IGwgMjUuMzg4MjY4LDQzLjUzNTE1NjMgbCAyNC42OTc0NDc3LDQzLjY4ODY3MTkgbCAyMy44NTMxMTE3LDQzLjY4ODY3MTkgbCAyMy4yMzkwNDkyLDQzLjkxODk0NTMgbCAyMi44NTUyNjAyLDQ0LjM3OTQ5MjIgbCAyMi4xNjQ0Mzk5LDQ0LjUzMzAwNzggbCAyMS41NTAzNzc0LDQ0Ljc2MzI4MTMgbCAyMS4xNjY1ODgzLDQ1LjIyMzgyODEgbCAyMC40NzU3NjgsNDUuMzc3MzQzOCBsIDE5LjYzMTQzMjEsNDUuMzc3MzQzOCBsIDE5LjAxNzM2OTYsNDUuNjA3NjE3MiBsIDE4LjYzMzU4MDUsNDYuMDY4MTY0MSBsIDE3Ljk0Mjc2MDIsNDYuMjIxNjc5NyBsIDE2LjI1NDA4ODMsNDYuMjIxNjc5NyBsIDE1LjY0MDAyNTgsNDYuNDUxOTUzMSBsIDE1LjI1NjIzNjcsNDYuOTEyNSBsIDE0Ljc5NTY4OTksNDYuOTEyNSBsIDE0LjQxMTkwMDgsNDYuNDUxOTUzMSBsIDEzLjcyMTA4MDUsNDYuMjIxNjc5NyBsIDkuNDk5NDAwOCw0Ni4yMjE2Nzk3IGwgOC44ODUzMzgzLDQ2LjA2ODE2NDEgbCA4LjUwMTU0OTI0LDQ1LjYwNzYxNzIgbCA3LjgxMDcyODkzLDQ1LjM3NzM0MzggbCA3LjE5NjY2NjQzLDQ1LjIyMzgyODEgbCA2LjgxMjg3NzM3LDQ0Ljc2MzI4MTMgbCA2LjEyMjA1NzA1LDQ0LjUzMzAwNzggbCA1LjUwNzk5NDU1LDQ0LjM3OTQ5MjIgbCA1LjEyNDIwNTQ5LDQzLjkxODk0NTMgbCA0LjY2MzY1ODYyLDQzLjUzNTE1NjMgbCA0LjI3OTg2OTU1LDQzLjA3NDYwOTQgbCAzLjgxOTMyMjY4LDQyLjY5MDgyMDMgbCAzLjQzNTUzMzYyLDQyLjIzMDI3MzQgbCAyLjk3NDk4Njc0LDQxLjg0NjQ4NDQgbCAyLjU5MTE5NzY4LDQxLjM4NTkzNzUgbCAyLjEzMDY1MDgsNDEuMDAyMTQ4NCBsIDEuOTc3MTM1MTgsNDAuMzg4MDg1OSBsIDEuNzQ2ODYxNzQsMzkuNjk3MjY1NiBsIDEuMjg2MzE0ODcsMzkuMzEzNDc2NiBsIDEuMTMyNzk5MjQsMzguNjk5NDE0MSBsIDEuMTMyNzk5MjQsMzcuODU1MDc4MSBsIDAuOTAyNTI1ODA0LDM3LjE2NDI1NzggbCAwLjQ0MTk3ODkyOSwzNi43ODA0Njg4IGwgMC4yODg0NjMzMDQsMzYuMTY2NDA2MyBsIDAuMjg4NDYzMzA0LDMzLjYzMzM5ODQgbCAwLjQ0MTk3ODkyOSwzMi45NDI1NzgxIGwgMC45MDI1MjU4MDQsMzIuNTU4Nzg5MSBsIDEuMTMyNzk5MjQsMzEuOTQ0NzI2NiBsIDEuMTMyNzk5MjQsMzEuMTAwMzkwNiBsIDEuMjg2MzE0ODcsMzAuNDA5NTcwMyBsIDEuNzQ2ODYxNzQsMzAuMDI1NzgxMyBsIDEuOTc3MTM1MTgsMjkuNDExNzE4OCBsIDIuMTMwNjUwOCwyOC43MjA4OTg0IGwgMi41OTExOTc2OCwyOC4zMzcxMDk0IGwgMi45NzQ5ODY3NCwyNy44NzY1NjI1IGwgMy40MzU1MzM2MiwyNy40OTI3NzM0IGwgMy44MTkzMjI2OCwyNy4wMzIyMjY2IGwgNC4yNzk4Njk1NSwyNi42NDg0Mzc1IGwgNC42NjM2NTg2MiwyNi4xODc4OTA2IGwgNS4yNzc3MjExMiwyNS45NTc2MTcyIGwgNi4xMjIwNTcwNSwyNS45NTc2MTcyIGwgNi44MTI4NzczNywyNS44MDQxMDE2IGwgNy4xOTY2NjY0MywyNS4zNDM1NTQ3IGwgNy44MTA3Mjg5MywyNS4xMTMyODEzIGwgOC42NTUwNjQ4NywyNS4xMTMyODEzIGwgOS4zNDU4ODUxOCwyNC45NTk3NjU2IGwgOS43Mjk2NzQyNCwyNC40OTkyMTg4IGwgMTAuMzQzNzM2NywyNC4yNjg5NDUzIGwgMTMuNzIxMDgwNSwyNC4yNjg5NDUzIGwgMTQuNDExOTAwOCwyNC4xMTU0Mjk3IGwgMTQuNzk1Njg5OSwyMy42NTQ4ODI4IGwgMTUuNDA5NzUyNCwyMy40MjQ2MDk0IGwgMTcuOTQyNzYwMiwyMy40MjQ2MDk0IGwgMTguNjMzNTgwNSwyMy4yNzEwOTM4IGwgMTkuMDE3MzY5NiwyMi44MTA1NDY5IGwgMTkuNjMxNDMyMSwyMi41ODAyNzM0IGwgMjMuODUzMTExNywyMi41ODAyNzM0IGwgMjQuNTQzOTMyMSwyMi40MjY3NTc4IGwgMjQuOTI3NzIxMSwyMS45NjYyMTA5IGwgMjUuMzg4MjY4LDIxLjU4MjQyMTkgbCAyNS43NzIwNTcxLDIxLjEyMTg3NSBsIDI2LjM4NjExOTYsMjAuODkxNjAxNiBsIDI3LjA3NjkzOTksMjAuNzM4MDg1OSBsIDI3LjMwNzIxMzMsMjAuMTI0MDIzNCBsIDI3LjMwNzIxMzMsMTguNDM1MzUxNiBsIDI3LjA3NjkzOTksMTcuNzQ0NTMxMyBsIDI2LjYxNjM5MywxNy4zNjA3NDIyIGwgMjYuNjE2MzkzLDE2LjkwMDE5NTMgbCAyNy4wNzY5Mzk5LDE2LjUxNjQwNjMgbCAyNy4wNzY5Mzk5LDE2LjA1NTg1OTQgbCAyNi42MTYzOTMsMTUuNjcyMDcwMyBsIDI2LjQ2Mjg3NzQsMTUuMDU4MDA3OCBsIDI2LjIzMjYwMzksMTQuMzY3MTg3NSBsIDI1Ljc3MjA1NzEsMTMuOTgzMzk4NCBsIDI1LjM4ODI2OCwxMy41MjI4NTE2IGwgMjQuOTI3NzIxMSwxMy4xMzkwNjI1IGwgMjQuNTQzOTMyMSwxMi42Nzg1MTU2IGwgMjMuODUzMTExNywxMi40NDgyNDIyIGwgMjMuMDA4Nzc1OCwxMi40NDgyNDIyIGwgMjIuMzk0NzEzMywxMi4yOTQ3MjY2IGwgMjIuMDEwOTI0MiwxMS44MzQxNzk3IGwgMjEuNTUwMzc3NCwxMS44MzQxNzk3IGwgMjEuMTY2NTg4MywxMi4yOTQ3MjY2IGwgMjAuNDc1NzY4LDEyLjQ0ODI0MjIgbCAxNC41NjU0MTY0LDEyLjQ0ODI0MjIgbCAxMy45NTEzNTM5LDEyLjY3ODUxNTYgbCAxMy41Njc1NjQ5LDEzLjEzOTA2MjUgbCAxMy4xMDcwMTgsMTMuNTIyODUxNiBsIDEyLjcyMzIyODksMTMuOTgzMzk4NCBsIDEyLjI2MjY4MjEsMTQuMzY3MTg3NSBsIDEyLjEwOTE2NjQsMTUuMDU4MDA3OCBsIDEyLjEwOTE2NjQsMTcuNTkxMDE1NiBsIDExLjg3ODg5MywxOC4yMDUwNzgxIGwgMTEuNDE4MzQ2MSwxOC41ODg4NjcyIGwgMTEuMDM0NTU3MSwxOS4wNDk0MTQxIGwgMTAuNTc0MDEwMiwxOS40MzMyMDMxIGwgMTAuMTkwMjIxMSwxOS44OTM3NSBsIDkuNDk5NDAwOCwyMC4wNDcyNjU2IGwgNy44MTA3Mjg5MywyMC4wNDcyNjU2IGwgNy4xOTY2NjY0MywxOS44OTM3NSBsIDYuODEyODc3MzcsMTkuNDMzMjAzMSBsIDYuMTIyMDU3MDUsMTkuMjAyOTI5NyBsIDUuNTA3OTk0NTUsMTkuMDQ5NDE0MSBsIDUuMTI0MjA1NDksMTguNTg4ODY3MiBsIDQuNjYzNjU4NjIsMTguMjA1MDc4MSBsIDQuNTEwMTQyOTksMTcuNTkxMDE1NiBsIDQuNTEwMTQyOTksMTYuNzQ2Njc5NyBsIDQuMjc5ODY5NTUsMTYuMDU1ODU5NCBsIDMuODE5MzIyNjgsMTUuNjcyMDcwMyBsIDMuNjY1ODA3MDUsMTUuMDU4MDA3OCBsIDMuNjY1ODA3MDUsMTQuMjEzNjcxOSBsIDMuODE5MzIyNjgsMTMuNTIyODUxNiBsIDQuMjc5ODY5NTUsMTMuMTM5MDYyNSBsIDQuNTEwMTQyOTksMTIuNTI1IGwgNC42NjM2NTg2MiwxMS44MzQxNzk3IGwgNS4xMjQyMDU0OSwxMS40NTAzOTA2IGwgNS41MDc5OTQ1NSwxMC45ODk4NDM4IGwgNS45Njg1NDE0MywxMC42MDYwNTQ3IGwgNi4zNTIzMzA0OSwxMC4xNDU1MDc4IGwgNi45NjYzOTI5OSw5LjkxNTIzNDM4IGwgNy42NTcyMTMzLDkuNzYxNzE4NzUgbCA4LjA0MTAwMjM3LDkuMzAxMTcxODggbCA4LjY1NTA2NDg3LDkuMDcwODk4NDQgbCA5LjM0NTg4NTE4LDguOTE3MzgyODEgbCA5LjcyOTY3NDI0LDguNDU2ODM1OTQgbCAxMC4zNDM3MzY3LDguMjI2NTYyNSBsIDEzLjcyMTA4MDUsOC4yMjY1NjI1IGwgMTQuNDExOTAwOCw4LjA3MzA0Njg4IGwgMTQuNzk1Njg5OSw3LjYxMjUgbCAxNS4yNTYyMzY3LDcuNjEyNSBsIDE1LjY0MDAyNTgsOC4wNzMwNDY4OCBsIDE2LjI1NDA4ODMsOC4yMjY1NjI1IGwgMTkuNjMxNDMyMSw4LjIyNjU2MjUgbCAyMC4zMjIyNTI0LDguNDU2ODM1OTQgbCAyMC43MDYwNDE0LDguOTE3MzgyODEgbCAyMS4zMjAxMDM5LDkuMDcwODk4NDQgbCAyMy44NTMxMTE3LDkuMDcwODk4NDQgbCAyNC41NDM5MzIxLDkuMzAxMTcxODggbCAyNC45Mjc3MjExLDkuNzYxNzE4NzUgbCAyNS41NDE3ODM2LDkuOTE1MjM0MzggbCAyNi4yMzI2MDM5LDEwLjE0NTUwNzggbCAyNi42MTYzOTMsMTAuNjA2MDU0NyBsIDI3LjA3NjkzOTksMTAuOTg5ODQzOCBsIDI3LjQ2MDcyODksMTEuNDUwMzkwNiBsIDI3LjkyMTI3NTgsMTEuODM0MTc5NyBsIDI4LjMwNTA2NDksMTIuMjk0NzI2NiBsIDI4Ljc2NTYxMTcsMTIuNjc4NTE1NiBsIDI5LjE0OTQwMDgsMTMuMTM5MDYyNSBsIDI5LjYwOTk0NzcsMTMuNTIyODUxNiBsIDI5Ljg0MDIyMTEsMTQuMjEzNjcxOSBsIDI5Ljk5MzczNjcsMTQuODI3NzM0NCBsIDMwLjQ1NDI4MzYsMTUuMjExNTIzNCBsIDMwLjY4NDU1NzEsMTUuOTAyMzQzOCBsIDMwLjY4NDU1NzEsMjAuMTI0MDIzNCBsIDMwLjgzODA3MjcsMjAuNzM4MDg1OSBsIDMxLjI5ODYxOTYsMjEuMTIxODc1IGwgMzEuNTI4ODkzLDIxLjgxMjY5NTMgbCAzMS41Mjg4OTMsMzUuMzIyMDcwMyBsIDMxLjY4MjQwODYsMzUuOTM2MTMyOCBsIDMyLjE0Mjk1NTUsMzYuMzE5OTIxOSBsIDMyLjUyNjc0NDYsMzYuNzgwNDY4OCBsIDMyLjk4NzI5MTQsMzcuMTY0MjU3OCBsIDMzLjM3MTA4MDUsMzcuNjI0ODA0NyBsIDMzLjk4NTE0MywzNy43NzgzMjAzIGwgMzQuNjc1OTYzMywzOC4wMDg1OTM4IGwgMzUuMDU5NzUyNCwzOC40NjkxNDA2IGwgMzUuNjczODE0OSwzOC42MjI2NTYzIGwgMzYuMzY0NjM1MiwzOC44NTI5Mjk3IGwgMzYuNzQ4NDI0MiwzOS4zMTM0NzY2IGwgMzcuMzYyNDg2NywzOS40NjY5OTIyIGwgMzguMDUzMzA3MSwzOS42OTcyNjU2IGwgMzguMjgzNTgwNSw0MC4zODgwODU5IGwgMzguMjgzNTgwNSw0Mi45MjEwOTM4IGwgMzguMDUzMzA3MSw0My41MzUxNTYzIGwgMzcuNTkyNzYwMiw0My45MTg5NDUzIHogbSA4Ny40MTgxODU3LDQwLjU0MTYwMTYgbCA4Ny4yNjQ2NzAxLDQxLjIzMjQyMTkgbCA4Ny4wMzQzOTY2LDQxLjg0NjQ4NDQgbCA4Ni41NzM4NDk3LDQyLjIzMDI3MzQgbCA4Ni4xOTAwNjA3LDQyLjY5MDgyMDMgbCA4NS43Mjk1MTM4LDQzLjA3NDYwOTQgbCA4NS4zNDU3MjQ3LDQzLjUzNTE1NjMgbCA4NC42NTQ5MDQ0LDQzLjY4ODY3MTkgbCA4NC4wNDA4NDE5LDQzLjkxODk0NTMgbCA4My42NTcwNTI5LDQ0LjM3OTQ5MjIgbCA4Mi45NjYyMzI2LDQ0LjUzMzAwNzggbCA3Mi44MzQyMDEzLDQ0LjUzMzAwNzggbCA3Mi4yMjAxMzg4LDQ0LjM3OTQ5MjIgbCA3MS44MzYzNDk3LDQzLjkxODk0NTMgbCA3MS4zNzU4MDI5LDQzLjUzNTE1NjMgbCA3MS4yMjIyODcyLDQyLjkyMTA5MzggbCA3MC45OTIwMTM4LDQyLjIzMDI3MzQgbCA3MC41MzE0NjY5LDQxLjg0NjQ4NDQgbCA3MC4zNzc5NTEzLDQxLjIzMjQyMTkgbCA3MC41MzE0NjY5LDQwLjU0MTYwMTYgbCA3MC45OTIwMTM4LDQwLjE1NzgxMjUgbCA3MS4zNzU4MDI5LDM5LjY5NzI2NTYgbCA3MS44MzYzNDk3LDM5LjMxMzQ3NjYgbCA3Mi4yMjAxMzg4LDM4Ljg1MjkyOTcgbCA3Mi42ODA2ODU3LDM4LjQ2OTE0MDYgbCA3My4wNjQ0NzQ3LDM4LjAwODU5MzggbCA3My42Nzg1MzcyLDM3Ljc3ODMyMDMgbCA3NC4zNjkzNTc2LDM3LjYyNDgwNDcgbCA3NC43NTMxNDY2LDM3LjE2NDI1NzggbCA3NS4yMTM2OTM1LDM2Ljc4MDQ2ODggbCA3NS41OTc0ODI2LDM2LjMxOTkyMTkgbCA3Ni4wNTgwMjk0LDM1LjkzNjEzMjggbCA3Ni4yODgzMDI5LDM1LjMyMjA3MDMgbCA3Ni4yODgzMDI5LDM0LjQ3NzczNDQgbCA3Ni40NDE4MTg1LDMzLjc4NjkxNDEgbCA3Ni45MDIzNjU0LDMzLjQwMzEyNSBsIDc3LjEzMjYzODgsMzIuNzg5MDYyNSBsIDc3LjEzMjYzODgsMTcuNTkxMDE1NiBsIDc3LjI4NjE1NDQsMTYuOTAwMTk1MyBsIDc3Ljc0NjcwMTMsMTYuNTE2NDA2MyBsIDc3Ljk3Njk3NDcsMTUuOTAyMzQzOCBsIDc3Ljc0NjcwMTMsMTUuMjExNTIzNCBsIDc3LjI4NjE1NDQsMTQuODI3NzM0NCBsIDc3LjEzMjYzODgsMTQuMjEzNjcxOSBsIDc2LjkwMjM2NTQsMTMuNTIyODUxNiBsIDc2LjQ0MTgxODUsMTMuMTM5MDYyNSBsIDc2LjA1ODAyOTQsMTIuNjc4NTE1NiBsIDc1LjM2NzIwOTEsMTIuNDQ4MjQyMiBsIDc0Ljc1MzE0NjYsMTIuMjk0NzI2NiBsIDc0LjM2OTM1NzYsMTEuODM0MTc5NyBsIDczLjY3ODUzNzIsMTEuNjAzOTA2MyBsIDczLjA2NDQ3NDcsMTEuNDUwMzkwNiBsIDcyLjY4MDY4NTcsMTAuOTg5ODQzOCBsIDcxLjk4OTg2NTQsMTAuNzU5NTcwMyBsIDY2LjA3OTUxMzgsMTAuNzU5NTcwMyBsIDY1LjQ2NTQ1MTMsMTAuOTg5ODQzOCBsIDY1LjA4MTY2MjIsMTEuNDUwMzkwNiBsIDY0LjM5MDg0MTksMTEuNjAzOTA2MyBsIDYzLjc3Njc3OTQsMTEuODM0MTc5NyBsIDYzLjYyMzI2MzgsMTIuNTI1IGwgNjMuMzkyOTkwNCwxMy4xMzkwNjI1IGwgNjIuOTMyNDQzNSwxMy41MjI4NTE2IGwgNjIuNTQ4NjU0NCwxMy45ODMzOTg0IGwgNjIuMDg4MTA3NiwxNC4zNjcxODc1IGwgNjEuNzA0MzE4NSwxNC44Mjc3MzQ0IGwgNjEuMjQzNzcxNiwxNS4yMTE1MjM0IGwgNjEuMDkwMjU2LDE1LjkwMjM0MzggbCA2MC44NTk5ODI2LDE2LjUxNjQwNjMgbCA2MC4zOTk0MzU3LDE2LjkwMDE5NTMgbCA2MC4wMTU2NDY2LDE3LjM2MDc0MjIgbCA1OS41NTUwOTk3LDE3Ljc0NDUzMTMgbCA1OS40MDE1ODQxLDE4LjQzNTM1MTYgbCA1OS4xNzEzMTA3LDE5LjA0OTQxNDEgbCA1OC43MTA3NjM4LDE5LjQzMzIwMzEgbCA1OC41NTcyNDgyLDIwLjEyNDAyMzQgbCA1OC4zMjY5NzQ3LDIwLjczODA4NTkgbCA1Ny44NjY0Mjc5LDIxLjEyMTg3NSBsIDU3LjcxMjkxMjIsMjEuODEyNjk1MyBsIDU3LjQ4MjYzODgsMjIuNDI2NzU3OCBsIDU3LjAyMjA5MTksMjIuODEwNTQ2OSBsIDU2Ljg2ODU3NjMsMjMuNTAxMzY3MiBsIDU2Ljg2ODU3NjMsMjguNTY3MzgyOCBsIDU3LjAyMjA5MTksMjkuMTgxNDQ1MyBsIDU3LjQ4MjYzODgsMjkuNTY1MjM0NCBsIDU3LjcxMjkxMjIsMzAuMjU2MDU0NyBsIDU3LjcxMjkxMjIsMzEuMTAwMzkwNiBsIDU3LjQ4MjYzODgsMzEuNzE0NDUzMSBsIDU3LjAyMjA5MTksMzIuMDk4MjQyMiBsIDU2Ljg2ODU3NjMsMzIuNzg5MDYyNSBsIDU3LjAyMjA5MTksMzMuNDAzMTI1IGwgNTcuNDgyNjM4OCwzMy43ODY5MTQxIGwgNTcuNzEyOTEyMiwzNC40Nzc3MzQ0IGwgNTcuNzEyOTEyMiwzNS4zMjIwNzAzIGwgNTcuODY2NDI3OSwzNS45MzYxMzI4IGwgNTguMzI2OTc0NywzNi4zMTk5MjE5IGwgNTguNzEwNzYzOCwzNi43ODA0Njg4IGwgNTkuMTcxMzEwNywzNy4xNjQyNTc4IGwgNTkuNTU1MDk5NywzNy42MjQ4MDQ3IGwgNjAuMTY5MTYyMiwzNy43NzgzMjAzIGwgNjAuODU5OTgyNiwzOC4wMDg1OTM4IGwgNjEuMjQzNzcxNiwzOC40NjkxNDA2IGwgNjEuNzA0MzE4NSwzOC44NTI5Mjk3IGwgNjIuMDg4MTA3NiwzOS4zMTM0NzY2IGwgNjIuNTQ4NjU0NCwzOS42OTcyNjU2IGwgNjIuNzc4OTI3OSw0MC4zODgwODU5IGwgNjIuOTMyNDQzNSw0MS4wMDIxNDg0IGwgNjMuMzkyOTkwNCw0MS4zODU5Mzc1IGwgNjMuMzkyOTkwNCw0MS44NDY0ODQ0IGwgNjIuOTMyNDQzNSw0Mi4yMzAyNzM0IGwgNjIuNzc4OTI3OSw0Mi45MjEwOTM4IGwgNjIuNTQ4NjU0NCw0My41MzUxNTYzIGwgNjEuODU3ODM0MSw0My42ODg2NzE5IGwgNjEuMjQzNzcxNiw0My45MTg5NDUzIGwgNjAuODU5OTgyNiw0NC4zNzk0OTIyIGwgNjAuMTY5MTYyMiw0NC41MzMwMDc4IGwgNTkuMzI0ODI2Myw0NC41MzMwMDc4IGwgNTguNzEwNzYzOCw0NC43NjMyODEzIGwgNTguMzI2OTc0Nyw0NS4yMjM4MjgxIGwgNTcuODY2NDI3OSw0NS4yMjM4MjgxIGwgNTcuNDgyNjM4OCw0NC43NjMyODEzIGwgNTcuMDIyMDkxOSw0NC43NjMyODEzIGwgNTYuNjM4MzAyOSw0NS4yMjM4MjgxIGwgNTUuOTQ3NDgyNiw0NS4zNzczNDM4IGwgNTUuMTAzMTQ2Niw0NS4zNzczNDM4IGwgNTQuNDg5MDg0MSw0NS4yMjM4MjgxIGwgNTQuMTA1Mjk1MSw0NC43NjMyODEzIGwgNTMuNjQ0NzQ4Miw0NC43NjMyODEzIGwgNTMuMjYwOTU5MSw0NS4yMjM4MjgxIGwgNTIuNTcwMTM4OCw0NS4zNzczNDM4IGwgNDkuMTkyNzk1MSw0NS4zNzczNDM4IGwgNDguNTc4NzMyNiw0NS4yMjM4MjgxIGwgNDguMTk0OTQzNSw0NC43NjMyODEzIGwgNDcuNTA0MTIzMiw0NC41MzMwMDc4IGwgNDYuODkwMDYwNyw0NC4zNzk0OTIyIGwgNDYuNTA2MjcxNiw0My45MTg5NDUzIGwgNDYuMDQ1NzI0Nyw0My41MzUxNTYzIGwgNDUuNjYxOTM1Nyw0My4wNzQ2MDk0IGwgNDUuMjAxMzg4OCw0Mi42OTA4MjAzIGwgNDUuMDQ3ODczMiw0Mi4wNzY3NTc4IGwgNDUuMDQ3ODczMiw0MS4yMzI0MjE5IGwgNDUuMjAxMzg4OCw0MC41NDE2MDE2IGwgNDUuNjYxOTM1Nyw0MC4xNTc4MTI1IGwgNDYuMDQ1NzI0NywzOS42OTcyNjU2IGwgNDYuNTA2MjcxNiwzOS4zMTM0NzY2IGwgNDYuODkwMDYwNywzOC44NTI5Mjk3IGwgNDcuNTA0MTIzMiwzOC42MjI2NTYzIGwgNDguMTk0OTQzNSwzOC40NjkxNDA2IGwgNDguNTc4NzMyNiwzOC4wMDg1OTM4IGwgNDkuMDM5Mjc5NCwzNy42MjQ4MDQ3IGwgNDkuNDIzMDY4NSwzNy4xNjQyNTc4IGwgNDkuODgzNjE1NCwzNi43ODA0Njg4IGwgNTAuMjY3NDA0NCwzNi4zMTk5MjE5IGwgNTAuNzI3OTUxMywzNS45MzYxMzI4IGwgNTAuOTU4MjI0NywzNS4zMjIwNzAzIGwgNTAuOTU4MjI0NywyNi4wMzQzNzUgbCA1MS4xMTE3NDA0LDI1LjM0MzU1NDcgbCA1MS41NzIyODcyLDI0Ljk1OTc2NTYgbCA1MS44MDI1NjA3LDI0LjM0NTcwMzEgbCA1MS44MDI1NjA3LDIzLjUwMTM2NzIgbCA1MS41NzIyODcyLDIyLjgxMDU0NjkgbCA1MS4xMTE3NDA0LDIyLjQyNjc1NzggbCA1MC45NTgyMjQ3LDIxLjgxMjY5NTMgbCA1MC45NTgyMjQ3LDE2Ljc0NjY3OTcgbCA1MC43Mjc5NTEzLDE2LjA1NTg1OTQgbCA1MC4yNjc0MDQ0LDE1LjY3MjA3MDMgbCA1MC4xMTM4ODg4LDE1LjA1ODAwNzggbCA0OS44ODM2MTU0LDE0LjM2NzE4NzUgbCA0OS40MjMwNjg1LDEzLjk4MzM5ODQgbCA0OS4wMzkyNzk0LDEzLjUyMjg1MTYgbCA0OC41Nzg3MzI2LDEzLjEzOTA2MjUgbCA0OC4xOTQ5NDM1LDEyLjY3ODUxNTYgbCA0Ny43MzQzOTY2LDEyLjI5NDcyNjYgbCA0Ny4zNTA2MDc2LDExLjgzNDE3OTcgbCA0Ni44OTAwNjA3LDExLjQ1MDM5MDYgbCA0Ni41MDYyNzE2LDEwLjk4OTg0MzggbCA0NS44MTU0NTEzLDEwLjc1OTU3MDMgbCA0NS4yMDEzODg4LDEwLjYwNjA1NDcgbCA0NS4yMDEzODg4LDEwLjE0NTUwNzggbCA0NS42NjE5MzU3LDkuNzYxNzE4NzUgbCA0NS44OTIyMDkxLDkuMTQ3NjU2MjUgbCA0NS44OTIyMDkxLDguMzAzMzIwMzEgbCA0Ni4wNDU3MjQ3LDcuNjEyNSBsIDQ2LjUwNjI3MTYsNy4yMjg3MTA5NCBsIDQ2Ljg5MDA2MDcsNi43NjgxNjQwNiBsIDQ3LjM1MDYwNzYsNi4zODQzNzUgbCA0Ny43MzQzOTY2LDUuOTIzODI4MTMgbCA0OC4zNDg0NTkxLDUuNjkzNTU0NjkgbCA0OS4wMzkyNzk0LDUuNTQwMDM5MDYgbCA0OS40MjMwNjg1LDUuMDc5NDkyMTkgbCA1MC4wMzcxMzEsNC44NDkyMTg3NSBsIDUwLjcyNzk1MTMsNS4wNzk0OTIxOSBsIDUxLjExMTc0MDQsNS41NDAwMzkwNiBsIDUxLjcyNTgwMjksNS42OTM1NTQ2OSBsIDUyLjU3MDEzODgsNS42OTM1NTQ2OSBsIDUzLjI2MDk1OTEsNS45MjM4MjgxMyBsIDUzLjY0NDc0ODIsNi4zODQzNzUgbCA1NC4yNTg4MTA3LDYuNTM3ODkwNjMgbCA1NC45NDk2MzEsNi43NjgxNjQwNiBsIDU1LjMzMzQyMDEsNy4yMjg3MTA5NCBsIDU1Ljc5Mzk2NjksNy42MTI1IGwgNTYuMTc3NzU2LDguMDczMDQ2ODggbCA1Ni42MzgzMDI5LDguNDU2ODM1OTQgbCA1Ny4wMjIwOTE5LDguOTE3MzgyODEgbCA1Ny40ODI2Mzg4LDkuMzAxMTcxODggbCA1Ny44NjY0Mjc5LDkuNzYxNzE4NzUgbCA1OC4zMjY5NzQ3LDEwLjE0NTUwNzggbCA1OC43MTA3NjM4LDEwLjYwNjA1NDcgbCA1OS4xNzEzMTA3LDEwLjYwNjA1NDcgbCA1OS41NTUwOTk3LDEwLjE0NTUwNzggbCA2MC4xNjkxNjIyLDkuOTE1MjM0MzggbCA2MS4wMTM0OTgyLDkuOTE1MjM0MzggbCA2MS43MDQzMTg1LDkuNzYxNzE4NzUgbCA2Mi4wODgxMDc2LDkuMzAxMTcxODggbCA2Mi41NDg2NTQ0LDguOTE3MzgyODEgbCA2Mi45MzI0NDM1LDguNDU2ODM1OTQgbCA2My4zOTI5OTA0LDguMDczMDQ2ODggbCA2My43NzY3Nzk0LDcuNjEyNSBsIDY0LjM5MDg0MTksNy4zODIyMjY1NiBsIDY1LjA4MTY2MjIsNy4yMjg3MTA5NCBsIDY1LjQ2NTQ1MTMsNi43NjgxNjQwNiBsIDY2LjA3OTUxMzgsNi41Mzc4OTA2MyBsIDY4LjYxMjUyMTYsNi41Mzc4OTA2MyBsIDY5LjMwMzM0MTksNi4zODQzNzUgbCA2OS42ODcxMzEsNS45MjM4MjgxMyBsIDcwLjMwMTE5MzUsNS42OTM1NTQ2OSBsIDczLjY3ODUzNzIsNS42OTM1NTQ2OSBsIDc0LjM2OTM1NzYsNS45MjM4MjgxMyBsIDc0Ljc1MzE0NjYsNi4zODQzNzUgbCA3NS4zNjcyMDkxLDYuNTM3ODkwNjMgbCA3Ni4wNTgwMjk0LDYuNzY4MTY0MDYgbCA3Ni40NDE4MTg1LDcuMjI4NzEwOTQgbCA3Ny4wNTU4ODEsNy4zODIyMjY1NiBsIDc3Ljc0NjcwMTMsNy42MTI1IGwgNzguMTMwNDkwNCw4LjA3MzA0Njg4IGwgNzguNTkxMDM3Miw4LjQ1NjgzNTk0IGwgNzguOTc0ODI2Myw4LjkxNzM4MjgxIGwgNzkuNDM1MzczMiw5LjMwMTE3MTg4IGwgNzkuODE5MTYyMiw5Ljc2MTcxODc1IGwgODAuMjc5NzA5MSwxMC4xNDU1MDc4IGwgODAuNTA5OTgyNiwxMC44MzYzMjgxIGwgODAuNjYzNDk4MiwxMS40NTAzOTA2IGwgODEuMTI0MDQ1MSwxMS44MzQxNzk3IGwgODEuMzU0MzE4NSwxMi41MjUgbCA4MS41MDc4MzQxLDEzLjEzOTA2MjUgbCA4MS45NjgzODEsMTMuNTIyODUxNiBsIDgyLjE5ODY1NDQsMTQuMjEzNjcxOSBsIDgyLjE5ODY1NDQsMTYuNzQ2Njc5NyBsIDgxLjk2ODM4MSwxNy4zNjA3NDIyIGwgODEuNTA3ODM0MSwxNy43NDQ1MzEzIGwgODEuNTA3ODM0MSwxOC4yMDUwNzgxIGwgODEuOTY4MzgxLDE4LjU4ODg2NzIgbCA4MS45NjgzODEsMTkuMDQ5NDE0MSBsIDgxLjUwNzgzNDEsMTkuNDMzMjAzMSBsIDgxLjM1NDMxODUsMjAuMTI0MDIzNCBsIDgxLjM1NDMxODUsMjEuODEyNjk1MyBsIDgxLjUwNzgzNDEsMjIuNDI2NzU3OCBsIDgxLjk2ODM4MSwyMi44MTA1NDY5IGwgODIuMTk4NjU0NCwyMy41MDEzNjcyIGwgODIuMTk4NjU0NCwyOC41NjczODI4IGwgODIuMzUyMTcwMSwyOS4xODE0NDUzIGwgODIuODEyNzE2OSwyOS41NjUyMzQ0IGwgODMuMDQyOTkwNCwzMC4yNTYwNTQ3IGwgODMuMDQyOTkwNCwzMS4xMDAzOTA2IGwgODIuODEyNzE2OSwzMS43MTQ0NTMxIGwgODIuMzUyMTcwMSwzMi4wOTgyNDIyIGwgODIuMTk4NjU0NCwzMi43ODkwNjI1IGwgODIuMzUyMTcwMSwzMy40MDMxMjUgbCA4Mi44MTI3MTY5LDMzLjc4NjkxNDEgbCA4My4wNDI5OTA0LDM0LjQ3NzczNDQgbCA4My4xOTY1MDYsMzUuMDkxNzk2OSBsIDgzLjY1NzA1MjksMzUuNDc1NTg1OSBsIDg0LjA0MDg0MTksMzUuOTM2MTMyOCBsIDg0LjY1NDkwNDQsMzYuMDg5NjQ4NCBsIDg1LjM0NTcyNDcsMzYuMzE5OTIxOSBsIDg1LjcyOTUxMzgsMzYuNzgwNDY4OCBsIDg2LjM0MzU3NjMsMzYuOTMzOTg0NCBsIDg3LjAzNDM5NjYsMzcuMTY0MjU3OCBsIDg3LjI2NDY3MDEsMzcuODU1MDc4MSBsIDg3LjQxODE4NTcsMzguNDY5MTQwNiBsIDg3Ljg3ODczMjYsMzguODUyOTI5NyBsIDg4LjEwOTAwNiwzOS41NDM3NSBsIDg3Ljg3ODczMjYsNDAuMTU3ODEyNSBsIDg3LjQxODE4NTcsNDAuNTQxNjAxNiB6IG0gMTMwLjQ4ODkyNCwzMi4wOTgyNDIyIGwgMTMwLjMzNTQwOCwzMi43ODkwNjI1IGwgMTMwLjMzNTQwOCwzNC40Nzc3MzQ0IGwgMTMwLjEwNTEzNSwzNS4wOTE3OTY5IGwgMTI5LjY0NDU4OCwzNS40NzU1ODU5IGwgMTI5LjQ5MTA3MiwzNi4xNjY0MDYzIGwgMTI5LjQ5MTA3MiwzNy44NTUwNzgxIGwgMTI5LjI2MDc5OSwzOC40NjkxNDA2IGwgMTI4LjgwMDI1MiwzOC44NTI5Mjk3IGwgMTI4LjY0NjczNiwzOS41NDM3NSBsIDEyOC40MTY0NjMsNDAuMTU3ODEyNSBsIDEyNy45NTU5MTYsNDAuNTQxNjAxNiBsIDEyNy41NzIxMjcsNDEuMDAyMTQ4NCBsIDEyNy4xMTE1OCw0MS4zODU5Mzc1IGwgMTI2LjcyNzc5MSw0MS44NDY0ODQ0IGwgMTI2LjI2NzI0NCw0Mi4yMzAyNzM0IGwgMTI1Ljg4MzQ1NSw0Mi42OTA4MjAzIGwgMTI1LjQyMjkwOCw0My4wNzQ2MDk0IGwgMTI1LjAzOTExOSw0My41MzUxNTYzIGwgMTI0LjU3ODU3Miw0My45MTg5NDUzIGwgMTI0LjE5NDc4Myw0NC4zNzk0OTIyIGwgMTIzLjUwMzk2Myw0NC41MzMwMDc4IGwgMTIyLjg4OTksNDQuNzYzMjgxMyBsIDEyMi41MDYxMTEsNDUuMjIzODI4MSBsIDEyMS44MTUyOTEsNDUuMzc3MzQzOCBsIDEyMC45NzA5NTUsNDUuMzc3MzQzOCBsIDEyMC4zNTY4OTIsNDUuNjA3NjE3MiBsIDExOS45NzMxMDMsNDYuMDY4MTY0MSBsIDExOS4yODIyODMsNDYuMjIxNjc5NyBsIDExNS4wNjA2MDMsNDYuMjIxNjc5NyBsIDExNC40NDY1NDEsNDYuMDY4MTY0MSBsIDExNC4wNjI3NTIsNDUuNjA3NjE3MiBsIDExMy4zNzE5MzIsNDUuMzc3MzQzOCBsIDExMi41Mjc1OTYsNDUuMzc3MzQzOCBsIDExMS45MTM1MzMsNDUuMjIzODI4MSBsIDExMS41Mjk3NDQsNDQuNzYzMjgxMyBsIDExMC44Mzg5MjQsNDQuNTMzMDA3OCBsIDExMC4yMjQ4NjEsNDQuMzc5NDkyMiBsIDEwOS44NDEwNzIsNDMuOTE4OTQ1MyBsIDEwOS4zODA1MjUsNDMuNTM1MTU2MyBsIDEwOC45OTY3MzYsNDMuMDc0NjA5NCBsIDEwOC41MzYxODksNDIuNjkwODIwMyBsIDEwOC4xNTI0LDQyLjIzMDI3MzQgbCAxMDcuNjkxODUzLDQxLjg0NjQ4NDQgbCAxMDcuMzA4MDY0LDQxLjM4NTkzNzUgbCAxMDYuODQ3NTE3LDQxLjAwMjE0ODQgbCAxMDYuNDYzNzI4LDQwLjU0MTYwMTYgbCAxMDYuMDAzMTgyLDQwLjE1NzgxMjUgbCAxMDUuNjE5MzkyLDM5LjY5NzI2NTYgbCAxMDUuMTU4ODQ2LDM5LjMxMzQ3NjYgbCAxMDUuMDA1MzMsMzguNjk5NDE0MSBsIDEwNC43NzUwNTcsMzguMDA4NTkzOCBsIDEwNC4zMTQ1MSwzNy42MjQ4MDQ3IGwgMTA0LjE2MDk5NCwzNy4wMTA3NDIyIGwgMTA0LjE2MDk5NCwzNi4xNjY0MDYzIGwgMTAzLjkzMDcyMSwzNS40NzU1ODU5IGwgMTAzLjQ3MDE3NCwzNS4wOTE3OTY5IGwgMTAzLjMxNjY1OCwzNC40Nzc3MzQ0IGwgMTAzLjMxNjY1OCwzMi43ODkwNjI1IGwgMTAzLjQ3MDE3NCwzMi4wOTgyNDIyIGwgMTAzLjkzMDcyMSwzMS43MTQ0NTMxIGwgMTA0LjE2MDk5NCwzMS4xMDAzOTA2IGwgMTA0LjE2MDk5NCwyNy43MjMwNDY5IGwgMTA0LjMxNDUxLDI3LjAzMjIyNjYgbCAxMDQuNzc1MDU3LDI2LjY0ODQzNzUgbCAxMDUuMDA1MzMsMjYuMDM0Mzc1IGwgMTA0Ljc3NTA1NywyNS4zNDM1NTQ3IGwgMTA0LjMxNDUxLDI0Ljk1OTc2NTYgbCAxMDQuMTYwOTk0LDI0LjM0NTcwMzEgbCAxMDQuMzE0NTEsMjMuNjU0ODgyOCBsIDEwNC43NzUwNTcsMjMuMjcxMDkzOCBsIDEwNS4wMDUzMywyMi42NTcwMzEzIGwgMTA1LjAwNTMzLDE1LjkwMjM0MzggbCAxMDQuNzc1MDU3LDE1LjIxMTUyMzQgbCAxMDQuMzE0NTEsMTQuODI3NzM0NCBsIDEwNC4xNjA5OTQsMTQuMjEzNjcxOSBsIDEwMy45MzA3MjEsMTMuNTIyODUxNiBsIDEwMy40NzAxNzQsMTMuMTM5MDYyNSBsIDEwMy4wODYzODUsMTIuNjc4NTE1NiBsIDEwMi4zOTU1NjQsMTIuNDQ4MjQyMiBsIDEwMS43ODE1MDIsMTIuMjk0NzI2NiBsIDEwMS4zOTc3MTMsMTEuODM0MTc5NyBsIDEwMC43MDY4OTIsMTEuNjAzOTA2MyBsIDk4LjE3Mzg4NDYsMTEuNjAzOTA2MyBsIDk3LjU1OTgyMjEsMTEuNDUwMzkwNiBsIDk3LjE3NjAzMzEsMTAuOTg5ODQzOCBsIDk2LjcxNTQ4NjIsMTAuNjA2MDU0NyBsIDk2LjMzMTY5NzEsMTAuMTQ1NTA3OCBsIDk1Ljg3MTE1MDMsOS43NjE3MTg3NSBsIDk1LjQ4NzM2MTIsOS4zMDExNzE4OCBsIDk1LjAyNjgxNDMsOC45MTczODI4MSBsIDk0Ljg3MzI5ODcsOC4zMDMzMjAzMSBsIDk1LjAyNjgxNDMsNy42MTI1IGwgOTUuNDg3MzYxMiw3LjIyODcxMDk0IGwgOTUuODcxMTUwMyw2Ljc2ODE2NDA2IGwgOTYuMzMxNjk3MSw2LjM4NDM3NSBsIDk2LjcxNTQ4NjIsNS45MjM4MjgxMyBsIDk3LjE3NjAzMzEsNS41NDAwMzkwNiBsIDk3LjU1OTgyMjEsNS4wNzk0OTIxOSBsIDk4LjE3Mzg4NDYsNC44NDkyMTg3NSBsIDEwMi4zOTU1NjQsNC44NDkyMTg3NSBsIDEwMy4wODYzODUsNS4wNzk0OTIxOSBsIDEwMy40NzAxNzQsNS41NDAwMzkwNiBsIDEwMy45MzA3MjEsNS41NDAwMzkwNiBsIDEwNC4zMTQ1MSw1LjA3OTQ5MjE5IGwgMTA0Ljc3NTA1Nyw0LjY5NTcwMzEzIGwgMTA1LjAwNTMzLDQuMDgxNjQwNjMgbCAxMDUuMTU4ODQ2LDMuMzkwODIwMzEgbCAxMDUuNjE5MzkyLDMuMDA3MDMxMjUgbCAxMDUuODQ5NjY2LDIuMzkyOTY4NzUgbCAxMDYuMDAzMTgyLDEuNzAyMTQ4NDQgbCAxMDYuNjE3MjQ0LDEuNDcxODc1IGwgMTA3LjMwODA2NCwxLjMxODM1OTM4IGwgMTA3LjUzODMzOCwwLjcwNDI5Njg3NSBsIDEwNy42OTE4NTMsMC4wMTM0NzY1NjI1IGwgMTA4LjE1MjQsMC4wMTM0NzY1NjI1IGwgMTA4LjM4MjY3NCwwLjcwNDI5Njg3NSBsIDEwOC41MzYxODksMS4zMTgzNTkzOCBsIDEwOC45OTY3MzYsMS43MDIxNDg0NCBsIDEwOS4zODA1MjUsMi4xNjI2OTUzMSBsIDEwOS44NDEwNzIsMi41NDY0ODQzOCBsIDExMC4wNzEzNDYsMy4yMzczMDQ2OSBsIDExMC4wNzEzNDYsNC4wODE2NDA2MyBsIDExMC4yMjQ4NjEsNC42OTU3MDMxMyBsIDExMC42ODU0MDgsNS4wNzk0OTIxOSBsIDExMS4wNjkxOTcsNS41NDAwMzkwNiBsIDExMS42ODMyNiw1LjY5MzU1NDY5IGwgMTEyLjM3NDA4LDUuOTIzODI4MTMgbCAxMTIuNzU3ODY5LDYuMzg0Mzc1IGwgMTEzLjIxODQxNiw2LjM4NDM3NSBsIDExMy42MDIyMDUsNS45MjM4MjgxMyBsIDExNC4yMTYyNjcsNS42OTM1NTQ2OSBsIDExNS45MDQ5MzksNS42OTM1NTQ2OSBsIDExNi41OTU3Niw1LjkyMzgyODEzIGwgMTE2Ljk3OTU0OSw2LjM4NDM3NSBsIDExNy41OTM2MTEsNi41Mzc4OTA2MyBsIDEyMi42NTk2MjcsNi41Mzc4OTA2MyBsIDEyMy4zNTA0NDcsNi43NjgxNjQwNiBsIDEyMy43MzQyMzYsNy4yMjg3MTA5NCBsIDEyNC4xOTQ3ODMsNy42MTI1IGwgMTI0LjU3ODU3Miw4LjA3MzA0Njg4IGwgMTI1LjAzOTExOSw4LjQ1NjgzNTk0IGwgMTI1LjI2OTM5Miw5LjE0NzY1NjI1IGwgMTI1LjAzOTExOSw5Ljc2MTcxODc1IGwgMTI0LjU3ODU3MiwxMC4xNDU1MDc4IGwgMTI0LjE5NDc4MywxMC42MDYwNTQ3IGwgMTIzLjczNDIzNiwxMC45ODk4NDM4IGwgMTIzLjM1MDQ0NywxMS40NTAzOTA2IGwgMTIyLjY1OTYyNywxMS42MDM5MDYzIGwgMTEzLjM3MTkzMiwxMS42MDM5MDYzIGwgMTEyLjc1Nzg2OSwxMS44MzQxNzk3IGwgMTEyLjM3NDA4LDEyLjI5NDcyNjYgbCAxMTEuOTEzNTMzLDEyLjY3ODUxNTYgbCAxMTEuNTI5NzQ0LDEzLjEzOTA2MjUgbCAxMTEuMDY5MTk3LDEzLjUyMjg1MTYgbCAxMTAuNjg1NDA4LDEzLjk4MzM5ODQgbCAxMTAuMjI0ODYxLDE0LjM2NzE4NzUgbCAxMTAuMDcxMzQ2LDE1LjA1ODAwNzggbCAxMTAuMDcxMzQ2LDE2Ljc0NjY3OTcgbCAxMTAuMjI0ODYxLDE3LjM2MDc0MjIgbCAxMTAuNjg1NDA4LDE3Ljc0NDUzMTMgbCAxMTAuNjg1NDA4LDE4LjIwNTA3ODEgbCAxMTAuMjI0ODYxLDE4LjU4ODg2NzIgbCAxMTAuMDcxMzQ2LDE5LjI3OTY4NzUgbCAxMTAuMDcxMzQ2LDI4LjU2NzM4MjggbCAxMTAuMjI0ODYxLDI5LjE4MTQ0NTMgbCAxMTAuNjg1NDA4LDI5LjU2NTIzNDQgbCAxMTAuOTE1NjgyLDMwLjI1NjA1NDcgbCAxMTAuOTE1NjgyLDMxLjEwMDM5MDYgbCAxMTEuMDY5MTk3LDMxLjcxNDQ1MzEgbCAxMTEuNTI5NzQ0LDMyLjA5ODI0MjIgbCAxMTEuNzYwMDE3LDMyLjc4OTA2MjUgbCAxMTEuNzYwMDE3LDMzLjYzMzM5ODQgbCAxMTEuOTEzNTMzLDM0LjI0NzQ2MDkgbCAxMTIuMzc0MDgsMzQuNjMxMjUgbCAxMTIuNjA0MzUzLDM1LjMyMjA3MDMgbCAxMTIuNzU3ODY5LDM1LjkzNjEzMjggbCAxMTMuMzcxOTMyLDM2LjA4OTY0ODQgbCAxMTQuMDYyNzUyLDM2LjMxOTkyMTkgbCAxMTQuNDQ2NTQxLDM2Ljc4MDQ2ODggbCAxMTUuMDYwNjAzLDM2LjkzMzk4NDQgbCAxMTguNDM3OTQ3LDM2LjkzMzk4NDQgbCAxMTkuMTI4NzY3LDM2Ljc4MDQ2ODggbCAxMTkuNTEyNTU3LDM2LjMxOTkyMTkgbCAxMTkuOTczMTAzLDM1LjkzNjEzMjggbCAxMjAuMzU2ODkyLDM1LjQ3NTU4NTkgbCAxMjAuODE3NDM5LDM1LjA5MTc5NjkgbCAxMjEuMDQ3NzEzLDM0LjQ3NzczNDQgbCAxMjEuMjAxMjI4LDMzLjc4NjkxNDEgbCAxMjEuNjYxNzc1LDMzLjQwMzEyNSBsIDEyMS44OTIwNDksMzIuNzg5MDYyNSBsIDEyMS44OTIwNDksMzEuOTQ0NzI2NiBsIDEyMi4wNDU1NjQsMzEuMjUzOTA2MyBsIDEyMi41MDYxMTEsMzAuODcwMTE3MiBsIDEyMi43MzYzODUsMzAuMjU2MDU0NyBsIDEyMi44ODk5LDI5LjU2NTIzNDQgbCAxMjMuMzUwNDQ3LDI5LjE4MTQ0NTMgbCAxMjMuNzM0MjM2LDI4LjcyMDg5ODQgbCAxMjQuMzQ4Mjk5LDI4LjQ5MDYyNSBsIDEyNS4wMzkxMTksMjguMzM3MTA5NCBsIDEyNS40MjI5MDgsMjcuODc2NTYyNSBsIDEyNi4wMzY5NzEsMjcuNjQ2Mjg5MSBsIDEyNi44ODEzMDcsMjcuNjQ2Mjg5MSBsIDEyNy41NzIxMjcsMjcuODc2NTYyNSBsIDEyNy45NTU5MTYsMjguMzM3MTA5NCBsIDEyOC41Njk5NzgsMjguNDkwNjI1IGwgMTI5LjI2MDc5OSwyOC43MjA4OTg0IGwgMTI5LjY0NDU4OCwyOS4xODE0NDUzIGwgMTMwLjEwNTEzNSwyOS41NjUyMzQ0IGwgMTMwLjMzNTQwOCwzMC4yNTYwNTQ3IGwgMTMwLjQ4ODkyNCwzMC44NzAxMTcyIGwgMTMwLjk0OTQ3MSwzMS4yNTM5MDYzIGwgMTMwLjk0OTQ3MSwzMS43MTQ0NTMxIGwgMTMwLjQ4ODkyNCwzMi4wOTgyNDIyIHpcIjtcbnZhciBsZXR0cmVTID0gXCJtIDEzOC44MTE1MTcsMzUuMTA2MTAxMSBjIDEzOC42Mzc4MSwzNS43MzY1MTExIDEzOC42MDY1NjQsMzUuODAxMDYxNCAxMzguMTk4NDU1LDM2LjIwOTE3MDQgYyAxMzguMDU2NzUzLDM2LjM1MDg3MjIgMTM3Ljk5OTE1NSwzNy4wMjE1OTI4IDEzNy45NDg3OTgsMzcuMjIzMDE5OSBjIDEzNy44Njk1MzcsMzcuNTQwMDYyMiAxMzcuOTI2NDExLDM3LjgyOTIzMzEgMTM3Ljk0ODc5OCwzOC4wNzc1NDc1IGMgMTM3Ljk4MzA3OCwzOC40NTc3ODA2IDEzNy45Mjc3NjYsMzguOTQxNTM2NyAxMzguMTI2MDM3LDM5LjIyMTc0NzkgYyAxMzguMjQwNzEzLDM5LjM4MzgxNjQgMTM4LjQzODU5NiwzOS40MzM4MzUgMTM4LjU2MDU0NCwzOS42MjAxOTY0IGMgMTM4Ljc0NzA2NiwzOS45MDUyNDA4IDEzOC43OTU5MSw0MC40NTg0NTI3IDEzOC45NTE1OTksNDAuODMzOTMzIGMgMTM5LjA2MTEzMSw0MS4wOTgwOTIyIDEzOS4yNjM5MTgsNDEuMjA0NTk1OCAxMzkuNDI5NTU2LDQxLjM3MDIzNDMgYyAxMzkuODYwODE5LDQxLjgwMTQ5NjUgMTM5LjM4NzcxMiw0My44OTY5NDMyIDEzOS42OTkzMjQsNDQuODk4OTgwNiBjIDEzOS43OTg1MTgsNDUuMjE3OTU1NyAxNDAuMDU4NTMzLDQ1LjMzNjg1MTUgMTQwLjI0OTM4Miw0NS41MTYyNDM4IGMgMTQwLjQ2ODM3Miw0NS43MjIwODc3IDE0MC42Njc2OTgsNDUuOTUyNzQwNSAxNDAuODgzNzE0LDQ2LjIwNjQxOTYgYyAxNDEuMjg3OTg3LDQ2LjY4MTE3OTEgMTQxLjM2MzU2Myw0Ni44MDQyNjYyIDE0MS43NTg3MzMsNDYuODg4OTQ1MiBjIDE0Mi4xNTM5MDMsNDYuOTczNjI0MSAxNDIuMzg5NTcxLDQ2Ljg4ODk0NTIgMTQyLjU3NzMsNDYuODg4OTQ1MiBjIDE0Mi45MTkwMTgsNDYuODg4OTQ1MiAxNDMuMDYyNTg2LDQ2LjU0NDc5MjQgMTQzLjE3MDA1NSw0Ni40MzczMjI5IGMgMTQzLjUyNDQsNDYuMDgyOTc3NiAxNDQuMzAyNjEsNDYuMDg0OTMyNyAxNDQuNDU2MTA3LDQ1Ljk1NzY5MjggYyAxNDQuNjA5NjAzLDQ1LjgzMDQ1MjggMTQ0LjUzODI1Miw0NS43ODg4OTggMTQ0LjY4MDE2OSw0NS42NDY5ODExIGMgMTQ1LjIzNTk0Niw0NS4wOTEyMDQxIDE0OC4zMTM2OCw0NS4yMTEzNDc4IDE0OS4yMjE2NjEsNDUuMzkyOTQ0MiBjIDE0OS43OTM1MjcsNDUuNTA3MzE3MiAxNDkuNTkyMjkzLDQ1Ljk0MDY3OTQgMTUwLjMxMTM0NCw0Ni4wODQ0ODk2IGMgMTUwLjgzNjk0NCw0Ni4xODk2MDk2IDE1MS43NjUwMDYsNDYuMDcwMzc5NSAxNTIuNDAwMTAxLDQ2LjIwNjQxOTYgYyAxNTMuMDM1MTk1LDQ2LjM0MjQ1OTYgMTUyLjg0MzI1OSw0Ni42OTg5NjQ5IDE1My4zMDM3OTcsNDYuODU1MzU2MyBjIDE1NC4yNDEzMDYsNDcuMTczNzE5IDE1NC45NDkwNjUsNDcuMDcyMDg4MSAxNTUuNDIwMzI5LDQ2LjgwNDI2NjIgYyAxNTUuNzAxMTUzLDQ2LjY0NDY3MjMgMTU1Ljg1NjU4Miw0Ni4zMDE4ODUyIDE1Ni4yMzg0NDQsNDYuMjA2NDE5NiBjIDE1Ny43OTU1MTksNDUuODE3MTUwOSAxNTkuNzQ1Mzk3LDQ2LjUwMTI2MyAxNjEuMzc2MTA3LDQ1Ljk1NzY5MjggYyAxNjEuNTU2MTY1LDQ1Ljg5NzY3MzQgMTYxLjU5NTUwMSw0NS42NTQyMTU5IDE2MS44NTU5NTcsNDUuNTE2MjQzOCBjIDE2Mi4yMTQzNTcsNDUuMzI2Mzg3MiAxNjIuNzYwNDUsNDUuMjA0NTcyMiAxNjIuOTk5MTI3LDQ1LjEyNTAxMzMgYyAxNjMuMzAxNDY5LDQ1LjAyNDIzMjYgMTYzLjI1NDIxMSw0NC44MTc2ODQzIDE2My41NDI3OTgsNDQuNjczMzkxIGMgMTYzLjgxMTkyMSw0NC41Mzg4MjkzIDE2NC4xOTg3NTEsNDQuNTA0MDI5NCAxNjQuNTM3NDY5LDQ0LjM2Mjg5OTEgYyAxNjQuODc2MTg3LDQ0LjIyMTc2ODcgMTY1LjI0NTYyOCw0My43OTU4NzE4IDE2NS41OTU5Niw0My40NDU1NDA3IGMgMTY1LjY5NDA3NCw0My4zNDc0MjYzIDE2Ni4wODk0MjcsNDIuOTUxNTc3IDE2Ni40MTQwMzEsNDIuNTg0NjMzNyBjIDE2Ni43Mzg2MzUsNDIuMjE3NjkwMyAxNjYuNjQ5Mzk4LDQxLjY5NTgzODEgMTY2LjgyMzgxLDQxLjM3MDIzNDMgYyAxNjYuOTQyNTg3LDQxLjE0ODQ5MjggMTY3LjIzMTcyOCw0MS4wNzE3MjAyIDE2Ny4zNDYxNjQsNDAuODMzOTMzIGMgMTY3LjUxNDc5Nyw0MC40ODM1Mjg2IDE2Ny40NzA1NjQsMzkuOTU4NDMzOSAxNjcuNjQyNTQyLDM5LjYyMDE5NjQgYyAxNjcuNzUyMzMyLDM5LjQwNDI2ODEgMTY4LjA2NDc1NiwzOS4zMzA0NjQ2IDE2OC4xNjQ3MzIsMzkuMTI2MjMyOCBjIDE2OC4zMzc1NzYsMzguNzczMTQzMiAxNjguMzgyMjM1LDM4LjIxOTI0MjUgMTY4LjUxNzM5NSwzNy45NTQ4Mzc1IGMgMTY4LjY2MzgwMiwzNy42Njg0Mjk0IDE2OC45NDYzNTMsMzcuNTM2Njk1MSAxNjkuMDExMzU4LDM3LjQwNjY4NTcgYyAxNjkuMzIyMjA3LDM2Ljc4NDk4ODEgMTY5LjEzODE5MSwzNS4yMjQxMTY2IDE2OS4zNTAwNzQsMzQuNTg4NDY5MSBjIDE2OS40Njk3NzQsMzQuMjI5MzY5MyAxNjkuODE1ODEsMzQuMjkwNjkwMiAxNjkuOTcxMDU4LDMzLjc2ODUwMjcgYyAxNzAuMTI2MzA2LDMzLjI0NjMxNTIgMTY5LjY3ODYxMywzMi43Njc2MzU0IDE2OS40MDkyOSwzMi40OTgzMTI5IGMgMTY5LjIwOTMwNSwzMi4yOTgzMjczIDE2OS4zNDEwMzUsMzEuNDQzNzY1MyAxNjguNTczODUsMzAuODQ3MDY1OSBjIDE2OC4zMTk4MTEsMzAuNjQ5NDgwNSAxNjguNDU2NTU1LDMwLjAwODY2NjcgMTY4LjE2NDczMiwyOS42MDUzMjM0IGMgMTY3Ljk5MzQzNywyOS4zNjg1NjkxIDE2Ny43MjQxNjksMjkuMjg1NjczOSAxNjcuNjQyNTQyLDI5LjA0MDc5NDYgYyAxNjcuNDU1ODk2LDI4LjQ4MDg1NjggMTY3LjYwOTI3NSwyOC40Mjc5NjIyIDE2Ny4zNDYxNjQsMjcuOTgyMzA0IGMgMTY3LjEyNDE5NCwyNy42MDYzMjkzIDE2Ni42MjI3OCwyNy4xNDU5NzU5IDE2Ni4yODc2NzMsMjYuODU1NDA4OCBjIDE2NS45NDQ2MzUsMjYuNTU3OTY0MyAxNjUuNzExMDU0LDI2LjI0ODQ5NDcgMTY1LjQ1NDMzMSwyNi4wMDY0NTMxIGMgMTY1LjE0NDYyOCwyNS43MTQ0NjEzIDE2NS4wNDUwNjEsMjUuNDg0NTk0IDE2NC42Nzg3NjYsMjUuMjcxMTAzNCBjIDE2NC4zNDIyNywyNS4wNzQ5ODA1IDE2My45MjIzOTgsMjUuMDQ2NjYxNyAxNjMuNTQyNzk4LDI0Ljg3NzM5NTUgYyAxNjMuMTk0NDc4LDI0LjcyMjA3NzMgMTYzLjA4Mzk3MiwyNC4zNTUyMDYyIDE2Mi41MDc1MzYsMjQuMjI4MTg3NyBjIDE2MS45MzExLDI0LjEwMTE2OTMgMTYxLjQ4NjUwNSwyNC4yMjc4OTM2IDE2MS4wMjM0NDIsMjQuMDcyOTQxOCBjIDE2MC42ODM5MjYsMjMuOTU5MzMxMyAxNjAuNjEzNzIyLDIzLjU2MDAxNjkgMTYwLjI0NzIxNCwyMy40Mzc4NDc4IGMgMTU5LjE0OTAzLDIzLjA3MTc4NjMgMTU2Ljk4NzEzMiwyMy43MTM0MjEgMTU1LjkyNzkwOCwyMy4xODM4MDkxIGMgMTU1LjcwNjY1NiwyMy4wNzMxODI5IDE1NS40OTEwNjEsMjIuNjI1NDMwOCAxNTQuOTEyNDE4LDIyLjU0ODcxMzMgYyAxNTQuMzMzNzc2LDIyLjQ3MTk5NTcgMTUzLjgxNTI2OSwyMi41NDA4ODk5IDE1My40Nzc3NDMsMjIuMzk1MzczNCBjIDE1My4xNjI4NTcsMjIuMjU5NjE3NyAxNTMuMDU5NzYzLDIxLjg2Njk2MTIgMTUyLjkwNzY4LDIxLjgyODk0MDMgYyAxNTAuOTU0NDQ0LDIxLjM0MDYzMTQgMTQ5LjU2MzQ2NiwyMS45NTUwOTI4IDE0OC4zNTg5NDIsMjEuNDM3MTQxMyBjIDE0Ny45NzU5ODcsMjEuMjcyNDY5IDE0Ny44OTg4MjYsMjEuMDk1MDUyNyAxNDcuMjA4NjUsMjAuNDE3NjE4MyBjIDE0Ni41MTg0NzUsMTkuNzQwMTgzOSAxNDUuOTg1NjY4LDE5LjMzMDkwMDIgMTQ1LjczODgyNiwxOC45NjM5NTY4IGMgMTQ1LjQ5MTk4NCwxOC41OTcwMTM0IDE0NS41ODM1ODEsMTguMTQ1MzkwMiAxNDUuMzY4MTgyLDE3Ljc5MjU1OTcgYyAxNDUuMTUyNzgzLDE3LjQzOTcyOTIgMTQ1LjE2NTkyNiwxNy43OTg0OTUzIDE0NC44OTYwOTQsMTYuOTg4OTk5NSBjIDE0NC44MTA1OTMsMTYuNzMyNDk1MiAxNDUuMjg3NTk5LDE2LjYwMTQ0MzQgMTQ1LjM2ODE4MiwxNi4yOTg4MjM3IGMgMTQ1LjYzOTM3MiwxNS4yODA0MDY3IDE0NS42NDM5NTcsMTUuMTg1MDk2NiAxNDYuMTY1OTgyLDE0LjY1OTQyNDYgYyAxNDYuMzYwMTI3LDE0LjQ2MzkyMjggMTQ2LjQ5MDYyOSwxMy40OTc2NzY1IDE0Ni43MjYwOSwxMy4yNjIyMTY0IGMgMTQ3LjE4NzI3MSwxMi44MDEwMzUzIDE0Ny43NzA0NjcsMTIuNDQyNzA1IDE0OC4wMjQ1MDYsMTEuNzgwMzI4MyBjIDE0OC4yNzg1NDUsMTEuMTE3OTUxNSAxNDcuOTgyMTY3LDExLjA3OTEwNTMgMTQ4LjM1ODk0MiwxMC45MDUzMDg5IGMgMTQ4LjczNTcxNywxMC43MzE1MTI1IDE0OC43NzEzMDksMTEuMTYwOTQyIDE0OS4yMjE2NjEsMTEuMzk5MjcxNiBjIDE0OS41MDQ1NjcsMTEuNTQ4OTg3MyAxNDkuODIyMzE0LDExLjUyMDM2ODcgMTUwLjE5Nzk0MSwxMS41NTQ1MTY2IGMgMTUxLjA3NDQ0MiwxMS42MzQxOTg1IDE1My40MzQwMzgsMTEuNDYwODM4NCAxNTMuNjEzMzQsMTEuNjQwMTQwOSBjIDE1My43NjY1MDgsMTEuNzkzMzA4NSAxNTMuOTcxNTc3LDEyLjA4NjQ2NTQgMTU0LjIwNjA5NiwxMi4yMDM3MjQ5IGMgMTU1LjI1NTgxLDEyLjcyODU4MTggMTU2LjgzMzk4OCwxMi4wMjg3MjYxIDE1Ny44NjE0MiwxMi41NDI0NDIxIGMgMTU3Ljk4NDcyOCwxMi42MDQwOTYyIDE1OC4yMzMzNzMsMTIuOTU4MzQ5MyAxNTguMzQxNDY4LDEzLjAxMjM5NjkgYyAxNTkuMDI0Mjg0LDEzLjM1MzgwNDggMTU5LjI3ODIwNCwxMy4yNjE1Nzg4IDE1OS42Mzk2ODUsMTMuNDU4ODc1OCBjIDE1OS45NDIzMDgsMTMuNjI0MDQ3MyAxNjAuMTEzNTU5LDEzLjg5NjUzNjQgMTYwLjQ4NjQ3OCwxNC4yNzEwOTI3IGMgMTYwLjY4NjgwMSwxNC40NzIyOTQ2IDE2MC44ODY2NTIsMTQuNjg1OTgxIDE2MS4xODQ3LDE0Ljk4NDAyODUgYyAxNjEuNjEyNDYsMTUuNDExNzg4NSAxNjIuODcxMTA0LDE2LjkwNDcwMzMgMTYzLjQzNTYzMywxNy4yNTc1MzM4IGMgMTY0LjAwMDE2MiwxNy42MTAzNjQzIDE2NC4wMjEyODIsMTcuMzI1MjA1IDE2NC4zMzkzODcsMTcuNDExNTAyNiBjIDE2NS4yMTU3NjksMTcuNjQ5MjUyOSAxNjYuMDMzMTE1LDE3LjEwNTc1NjQgMTY2LjUwOTY5NywxNi4yOTg4MjM3IGMgMTY2Ljg2OTA3OCwxNS42OTAzMzE4IDE2Ni42OTMyMDksMTQuNTIzNzk5MSAxNjYuNjgyMTgyLDEzLjI2MjIxNjQgYyAxNjYuNjc0NzMzLDEyLjQxMDEyMzUgMTY2Ljc3NTU0OCwxMS4zMjMyNzg0IDE2Ni40MTQwMzEsMTAuOTYxNzYxNyBjIDE2Ni4xMjgwMjgsMTAuNjc1NzU5MiAxNjYuMTU5MDc2LDEwLjc5MDU3MDkgMTY1Ljk2MjQwNywxMC4zOTcyMzI5IGMgMTY1Ljg1NjIxMiwxMC4xODQ4NDQ0IDE2NS45MjQ4NzcsOS41OTI0MjIwNiAxNjUuNDU0MzMxLDkuMTY5MzgyOTcgYyAxNjQuOTgzNzg1LDguNzQ2MzQzODggMTY0Ljk3NDQ4MSw4LjM2MTgwNzk3IDE2NC40NTIyOTIsOC4yNjYxMzY4NSBjIDE2Mi41NzAyNjgsNy45MjEzMjgyNyAxNjAuODA0NzUxLDguMTU2NTA1MTMgMTU4Ljg0OTM0NSw4LjE1MzIzMTE2IGMgMTU4LjMwNDcxLDguMTUyMzE5MjcgMTU3LjUyMzUsOC4xOTc5MjQ3NSAxNTcuMDI4NzQsOC4wNzQyMzQ4OCBjIDE1Ni40Nzg2NzMsNy45MzY3MTc5NCAxNTYuNDEwNzU3LDcuNTM3OTUzNTIgMTU1LjkyNzkwOCw3LjM3NzAwNDEgYyAxNTUuNDQ4NjI0LDcuMjE3MjQyNTUgMTU0LjY1OTkwOSw3LjM3NzAwNDEgMTUzLjcwNzgsNy4yNzg4NzMxMyBjIDE1My4yNDExNjksNy4yMzA3Nzg4NSAxNTMuMTc3NDI0LDYuOTEzNjUxOTQgMTUyLjkwNzY4LDYuNzI4NDU3NiBjIDE1Mi42Mzc2ODgsNi41NDMwOTMxNCAxNTIuMjMwNTUsNi40MTc5NDA4NiAxNTEuNTY2OTI0LDYuNDc0NDE5NjYgYyAxNTEuMTIwNzgsNi41MTIzODkzMyAxNTAuODMyMDcsNi41NzMyMTIxOCAxNTAuMzcxOTUzLDcuMTIzNjI3NzEgYyAxNDkuODAyNzcxLDcuMzYzNTUyNDEgMTQ5LjI3OTY5Miw3LjM3MTA1NDMgMTQ5LjAxMjQzMSw3LjUwNDY4NDYxIGMgMTQ4LjcxNjA1NCw4LjQ5MjYwOTkxIDE0Ni44ODMyNiw3Ljk5NzAyMyAxNDYuNTE4NDc1LDguMzYxODA3OTYgYyAxNDYuMzk0ODIxLDguNDg1NDYxNjkgMTQ2LjMxNDk1Niw4LjY1MzgyMjI2IDE0Ni4xNjU5ODIsOC44MDI3OTY1OSBjIDE0NS45NDg2NzgsOS4wMjAxMDA4NyAxNDUuNDQwMTkxLDguODkwNzkxNTggMTQ1LjEzMDMwOCw5LjA0NTczMjc2IGMgMTQ0Ljg4MDQzLDkuMTcwNjcxOTQgMTQ0LjY2MDI2LDkuNDQ4NjM3NDYgMTQ0LjQ1NjEwNyw5LjY1Mjc5MDcyIGMgMTQ0LjM0Njg0Myw5Ljc2MjA1NDEyIDE0My43OTgxNzIsOS42NzA4NDQ5MiAxNDMuMzA0MzM1LDEwLjAyOTczMTYgYyAxNDIuODEwNDk4LDEwLjM4ODYxODQgMTQyLjE1Nzg2NSwxMS4wNDMyNTEyIDE0Mi4wODk4NSwxMS4xMzU1NTgyIGMgMTQxLjk1NzYyOCwxMS4zMTUwMDExIDE0MS42NDM0MjQsMTEuNTEwODc0MyAxNDEuNTQzOSwxMS43MDk5MjI0IGMgMTQxLjEwMjc1MiwxMi41OTIyMTc2IDE0MS4zNTc4MjIsMTIuNDcwMzgwOSAxNDEuMjM2NDkzLDEyLjY5ODAxNDQgYyAxNDEuMTc2NDA1LDEyLjgxMDc1MDUgMTQxLjA2MTM1OCwxMy4wNDc1Mzk2IDE0MC42NTAwMjIsMTMuNDU4ODc1OSBjIDE0MC42NDI0NTQsMTMuNDY2NDQzOSAxNDAuNDc1MzI2LDEzLjg2NzE0NDUgMTQwLjQ3OTQ0MSwxNC4yNzEwOTI0IGMgMTQwLjQ4MzIzMiwxNC42NDMzMzk1IDE0MC40OTYyNDQsMTUuMDA1NTgzMiAxNDAuNDA5NDI0LDE1LjE3OTIyMzUgYyAxNDAuMTE1NDAxLDE1Ljc2NzI2ODYgMTQwLjA3OTI5MiwxNS41Mjg4OTk3IDEzOS44NDgwMDcsMTUuOTkzNTM0NCBjIDEzOS43NTc5NDEsMTYuMTc0NDcwMSAxMzkuNTgwMzYsMTYuNjg4OTI3OCAxMzkuNjk5MzI0LDE2Ljk4ODk5OTIgYyAxMzkuODg1ODU0LDE3LjQ1OTQ5NzkgMTQwLjM3NDQ3NiwxNy45MDI4ODg2IDE0MC40MDk0MjQsMTguMDc3NjI3NCBjIDE0MC40NzY2OCwxOC40MTM5MDggMTQwLjQ3OTMzNCwxOC42NjI4MzMgMTQwLjQ3OTQ0MSwxOS4xMTcwMzk4IGMgMTQwLjQ3OTY4OCwyMC4xNzU4ODM0IDE0MS40NDI2LDIwLjI5MTc1MTYgMTQxLjkwNjA3MywyMS4wNTkzNzkyIGMgMTQyLjM2OTU0NywyMS44MjcwMDY3IDE0Mi4yMDU4ODQsMjIuMzMxMTk4NiAxNDIuMzY5NTQ3LDIyLjM5NTM3MzQgYyAxNDIuNzYwNjAzLDIyLjU0ODcxMzMgMTQzLjMyODg2NiwyMi40OTY2NTU2IDE0My42MDA2NDksMjIuNzY4NDM4IGMgMTQzLjg4MzgyMSwyMy4wNTE2MSAxNDMuODMyMzg1LDIzLjk5OTU0MDMgMTQ0LjA5MzA4OSwyNC4wNzI5NDE4IGMgMTQ0Ljc4NjczLDI0LjI2ODIzNjkgMTQ0LjkxODQ1LDI0LjMzMjY2MDQgMTQ1LjI2NjQyNiwyNC40NDg1Mjk1IGMgMTQ1LjYxNDQwMywyNC41NjQzOTg3IDE0NS42ODgwMjUsMjQuOTI2NDg3OSAxNDYuMTY1OTgyLDI1LjMzMjAyNjIgYyAxNDYuNjQzOTM4LDI1LjczNzU2NDYgMTQ2Ljg4ODU4NCwyNi4wMDY0NTMxIDE0Ny4yMDg2NSwyNi4zNzQ4NDEyIGMgMTQ3LjY4MTgyMiwyNi45MTk0NDk4IDE0OC40NDM1NTcsMjYuNTYzMTI4NSAxNDguOTQ1MjQ3LDI3LjMwMTc4ODggYyAxNDkuMzUzNTI5LDI3LjkwMjkxODUgMTUxLjg4OTY2NCwyNy4zMTExNTc4IDE1Mi43MTA5NywyNy43MjE4MTA4IGMgMTUzLjExOTc3NiwyNy45MjYyMTM4IDE1Mi45NDc4NiwyNy45MDQ1MjQ1IDE1My4xOTAxNjksMjguMTQ2ODMzNCBjIDE1My43OTA3MjEsMjguNzQ3Mzg1NiAxNTUuMzMyMzM3LDI4LjE4NTIwNyAxNTYuMDg1NjM2LDI4LjU2MTg1NjYgYyAxNTYuNDM3OTY2LDI4LjczODAyMTEgMTU2LjM2OTgwNSwyOC45Mjg0MzYzIDE1Ni45NDAxNjYsMjkuMjEzNjE2OCBjIDE1Ny4xNjQ5MjcsMjkuMzI1OTk3NCAxNTguMjAwMjM0LDI5LjIxMzYxNjggMTU4LjYzNDc0MSwyOS40MzA4Njk2IGMgMTU5LjA2OTI0OCwyOS42NDgxMjI0IDE1OS4wMzI1NSwyOS44NTk4ODk3IDE1OS4xOTU3NTYsMjkuOTE0MjkxOSBjIDE1OS41Nzk3OTcsMzAuMDQyMzA1NCAxNjAuMzg2NjE2LDMwLjE2Njk2NTQgMTYwLjYwOTU1NCwzMC4zODk5MDMyIGMgMTYwLjk5ODk2LDMwLjc3OTMwOTUgMTYxLjM1OTE4OCwzMS4yNTA4NjI4IDE2MS42NDczMTksMzEuNTE2NTAxMyBjIDE2Mi4xODkxMTQsMzIuMDE2MDAyNyAxNjIuNDM4MTYyLDMyLjMwNzIzNDkgMTYyLjQzODE2MiwzMy42MTA3MjA3IGMgMTYyLjQzODE2MiwzNC4wMTg2NjE5IDE2Mi40NTI4MzcsMzQuNzM3NDA4NSAxNjIuMjI2NjYxLDM0Ljk2MzU4NDUgYyAxNjIuMDk5NDE1LDM1LjA5MDgyOTYgMTYxLjM0MzA0LDM1LjE3NDA4NDggMTYxLjE4NDY5OSwzNS4yNTMyNTU1IGMgMTYwLjgwMzE1MywzNS40NDQwMjg3IDE2MC40MzY2MzYsMzUuODk5MDgxNSAxNjAuMTI2NTQ3LDM2LjIwOTE3MDQgYyAxNTkuODQxNzk0LDM2LjQ5MzkyMzMgMTU5LjE3OTE3NywzNy4yNDQyMjUxIDE1OC44NDkzNDQsMzcuNTI3MTcyNyBjIDE1OC41NjIzMjMsMzcuNzczMzk0NyAxNTcuOTQ1ODMsMzcuNjg4NzgzMSAxNTcuNjQ5ODU5LDM3Ljg0NTgxMjkgYyAxNTcuMjg4NTU5LDM4LjAzNzUwMzUgMTU3LjMwNDMzNCwzOC4yMjk0MjMyIDE1Ny4wMjg3MzksMzguMzY3MjIwMyBjIDE1Ni4yODgxMzUsMzguNzM3NTIyNyAxNTMuOTQ1NjI4LDM4LjM2NzIyMDMgMTUzLjQ3Nzc0MiwzOC43MDAzNDA1IGMgMTUyLjkwNzE4MywzOS4xMDY1NjA2IDE1My4wMjM1NDgsMzkuNDk2OTM3MSAxNTEuMjMzNjQ4LDM5LjMyMzEzMzQgYyAxNTAuNDE4NzQxLDM5LjI0NDAwMzkgMTUwLjU5NjEwNiwzOC44NTk2NjA1IDE1MC4xOTc5NDEsMzguNzAwMzQwNSBjIDE0OS43OTk3NzcsMzguNTQxMDIwNCAxNDguNDAzNDY4LDM4LjU0OTYxODMgMTQ4LjIyMTA3LDM4LjM2NzIyMDMgYyAxNDcuOTYyODE0LDM4LjEwODk2NDYgMTQ4LjA0NTg5OCwzOC4wNzY3Mjk5IDE0Ny43ODY1NjMsMzcuOTU0ODM3NSBjIDE0Ny40NzU1MjcsMzcuODA4NjQ0MSAxNDYuOTc3OTQ1LDM3LjcxMjU1MTcgMTQ2LjYwNjcyOSwzNy40Njc3MDc0IGMgMTQ1LjkyNjAwMiwzNy4wMTg3MTgxIDE0NS4yOTUzOTMsMzYuMTIyMjY4NiAxNDQuNTcxMjE2LDM1LjUyODQ0MjkgYyAxNDQuMzMwMzg3LDM1LjMzMDk2MjcgMTQzLjczMTE3LDM1LjE2OTA5MjkgMTQzLjMwNDMzNSwzNS4wMzYwMDI3IGMgMTQyLjg3NzUsMzQuOTAyOTEyNSAxNDMuMDgwOSwzNC4wNTc3Nzg0IDE0Mi44NzM5MzYsMzMuODE1NDI3MSBjIDE0Mi40NDUxNTQsMzMuMzEzMzMwMiAxNDEuNTc1NTcxLDMzLjMxMjQ2MDIgMTM5Ljk5NDQxNSwzMy42MTA3MjA3IGMgMTM5LjYyNDg3LDMzLjY4MDQyOTYgMTM5LjYxNzg0MywzNC4wOTQ1NzE0IDEzOC45NTE1OTksMzQuNTg4NDY5MSBsIDEzOC44MTE1MTcsMzUuMTA2MTAxMSB6XCI7XG4vLyB2YXIgcGF0aCA9IHBvaW50MSArIHBvaW50MiArIHBvaW50MyArIHBvaW50NCArIHBvaW50NSArIHBvaW50NiArIHBvaW50NyArIHBvaW50OCArIHBvaW50OTtcbnZhciBwYXRoID0gbGV0dHJlc0FOVCArIGxldHRyZVM7XG5cbm1vZHVsZS5leHBvcnRzID0gcGF0aDsiLCIndXNlIHN0cmljdCc7XG5cbnZhciBzcXJ0ID0gTWF0aC5zcXJ0O1xudmFyIHBvdyA9IE1hdGgucG93O1xuXG5mdW5jdGlvbiBzaWduKHgpIHtcblx0cmV0dXJuIHggPyB4IDwgMCA/IC0xIDogMSA6IDA7XG59XG5cbmZ1bmN0aW9uIHJhbmdlKHN0YXJ0LCBjb3VudCkge1xuICAgIHJldHVybiBBcnJheS5hcHBseSgwLCBBcnJheShjb3VudCkpLm1hcChmdW5jdGlvbiAoZWxlbWVudCwgaW5kZXgpIHtcbiAgICBcdHJldHVybiBpbmRleCArIHN0YXJ0XG4gICAgfSk7XG59XG5cbmZ1bmN0aW9uIGRpc3RhbmNlKGEsIGIpe1xuXHRyZXR1cm4gc3FydChwb3coYS54IC0gYi54LCAyKSArIHBvdyhhLnkgLSBiLnksIDIpKTtcbn1cblxuZnVuY3Rpb24gbm9ybSh2KXtcblx0cmV0dXJuIHNxcnQocG93KHYueCwgMikgKyBwb3codi55LCAyKSk7XG59XG5cbm1vZHVsZS5leHBvcnRzID0ge1xuXHRzaWduOiBzaWduLFxuXHRyYW5nZTogcmFuZ2UsXG5cdGRpc3RhbmNlOiBkaXN0YW5jZSxcblx0bm9ybTogbm9ybVxufSIsIid1c2Ugc3RyaWN0J1xuXG5mdW5jdGlvbiBWZWN0b3IoeCwgeSkge1xuICAgIHRoaXMueCA9IHg7ICAgICAgICAgICAgICAgIFxuICAgIHRoaXMueSA9IHk7XG59XG5cblZlY3Rvci5wcm90b3R5cGUubm9ybSA9IGZ1bmN0aW9uKCl7XG5cdHJldHVybiBNYXRoLnNxcnQodGhpcy54ICogdGhpcy54ICsgdGhpcy55ICogdGhpcy55KTtcbn1cblxuVmVjdG9yLnByb3RvdHlwZS5ub3JtYWxpemUgPSBmdW5jdGlvbigpe1xuXHR2YXIgbm9ybSA9IHRoaXMubm9ybSgpO1xuXHR0aGlzLnggPSB0aGlzLnggLyBub3JtO1xuXHR0aGlzLnkgPSB0aGlzLnkgLyBub3JtO1xufVxuXG5cblxubW9kdWxlLmV4cG9ydHMgPSBWZWN0b3I7IiwiJ3VzZSBzdHJpY3QnO1xuXG52YXIgX2FudENvbG9ueSA9IHJlcXVpcmUoJy4vaW5kZXguanMnKTtcblxudmFyIGNvbnRhaW5lciA9IGRvY3VtZW50LnF1ZXJ5U2VsZWN0b3IoJy5jb2xvbnknKTtcblxudmFyIG9wdGlvbnMgPSB7XG5cdHZlbG9jaXR5OiAwLjAwMSxcblx0bmJBbnRzOiA0MDAwLFxuXHR3ZWlnaHQ6IDEwLFxuXHRyZXBTaXplOiAwLjA1LFxuXHRyZXBTcGVlZDogMC4wMDIsXG5cdG5iU3RhcnQ6IDMwMCxcblx0bmJSYW5kOiAzMDBcblx0Ly8gb2JqIHBhciBkZWZhdXRcbn07XG5cbnZhciBhbnRDb2xvbnkgPSBfYW50Q29sb255KGNvbnRhaW5lciwgb3B0aW9ucyk7XG5cbndpbmRvdy5hZGRFdmVudExpc3RlbmVyKCdjbGljaycsIGZ1bmN0aW9uICgpe1xuXHQvLyBvcHRpb25zLnZlbG9jaXR5ID0gMC4wMDM7XG5cdG9wdGlvbnMubmJBbnRzID0gMjAwMDA7XG5cdC8vIG9wdGlvbnMud2VpZ2h0ID0gMTAwMDAwMDA7XG5cdC8vIG9wdGlvbnMucmVwU3BlZWQgPSAwLjAxO1xuXHQvLyBvcHRpb25zLnJlcFNpemUgPSAwLjE7XG5cblx0Ly8gYW50Q29sb255LmNoYW5nZU9wdGlvbnMob3B0aW9ucyk7XG5cdGFudENvbG9ueS5jaGFuZ2VPcHRpb25zKG9wdGlvbnMpO1xufSk7XG5cbiJdfQ==
