(function e(t,n,r){function s(o,u){if(!n[o]){if(!t[o]){var a=typeof require=="function"&&require;if(!u&&a)return a(o,!0);if(i)return i(o,!0);var f=new Error("Cannot find module '"+o+"'");throw f.code="MODULE_NOT_FOUND",f}var l=n[o]={exports:{}};t[o][0].call(l.exports,function(e){var n=t[o][1][e];return s(n?n:e)},l,l.exports,e,t,n,r)}return n[o].exports}var i=typeof require=="function"&&require;for(var o=0;o<r.length;o++)s(r[o]);return s})({"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/index.js":[function(require,module,exports){
'use strict';

var initRendering = require('./src/rendering.js');
var initializePoints = require('./src/initializePoints.js');
var createEdges = require('./src/createEdges.js');
var initAnts = require('./src/initializeAnts');

module.exports = function init(containerElement, options){

	var render, pointsInfos, edges, population, initVar;


	function _init(containerElement, options){
		pointsInfos = initializePoints(options.nbStart, options.nbRand);
		edges = createEdges(pointsInfos.points);
		population = initAnts(containerElement, pointsInfos, options);
		initVar = {
			pointsInfos: pointsInfos,
			edges: edges,
			population: population
		};
		render = initRendering(containerElement, initVar);
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
},{"./src/createEdges.js":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/createEdges.js","./src/initializeAnts":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/initializeAnts.js","./src/initializePoints.js":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/initializePoints.js","./src/rendering.js":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/rendering.js"}],"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/ich.js":[function(require,module,exports){
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
},{"robust-orientation":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/orientation.js","simplicial-complex":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/simplicial-complex/topology.js"}],"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-scale/node_modules/two-sum/two-sum.js":[function(require,module,exports){
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
},{}],"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-scale/robust-scale.js":[function(require,module,exports){
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
},{"two-product":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/two-product/two-product.js","two-sum":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-scale/node_modules/two-sum/two-sum.js"}],"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-subtract/robust-diff.js":[function(require,module,exports){
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
},{}],"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-sum/robust-sum.js":[function(require,module,exports){
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
},{}],"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/two-product/two-product.js":[function(require,module,exports){
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
},{}],"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/orientation.js":[function(require,module,exports){
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
},{"robust-scale":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-scale/robust-scale.js","robust-subtract":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-subtract/robust-diff.js","robust-sum":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-sum/robust-sum.js","two-product":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/two-product/two-product.js"}],"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/simplicial-complex/node_modules/bit-twiddle/twiddle.js":[function(require,module,exports){
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


},{}],"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/simplicial-complex/node_modules/union-find/index.js":[function(require,module,exports){
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
},{}],"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/simplicial-complex/topology.js":[function(require,module,exports){
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

},{"bit-twiddle":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/simplicial-complex/node_modules/bit-twiddle/twiddle.js","union-find":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/simplicial-complex/node_modules/union-find/index.js"}],"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/node_modules/uniq/uniq.js":[function(require,module,exports){
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

},{}],"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/triangulate.js":[function(require,module,exports){
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
},{"incremental-convex-hull":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/ich.js","uniq":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/node_modules/uniq/uniq.js"}],"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/parse-svg-path/index.js":[function(require,module,exports){

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

},{}],"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/ant.js":[function(require,module,exports){
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

            possibleEdges.splice(possibleEdges.indexOf(this.edge),1);

            // flip a coin and either take the smelliest path or a random one
            if (random() > 0.5){
                var smells = possibleEdges.map(function(e){return e.pheromon});
                var index = smells.indexOf(Math.max.apply(Math, smells));
                this.edge = possibleEdges[index];
            } 
            else
                this.edge = possibleEdges[floor(random()*possibleEdges.length)];

            // set the destination point, being edge.pt1 or edge.pt2
            this.destination = (this.origin == this.edge.pt1) ? this.edge.pt2 : this.edge.pt1;

            this.direction.x = this.destination.x - this.origin.x; 
            this.direction.y = this.destination.y - this.origin.y;

            this.direction.normalize();
        },

        move: function(){
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


},{"./mouse.js":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/mouse.js","./utilities.js":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/utilities.js","./vector.js":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/vector.js"}],"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/createEdges.js":[function(require,module,exports){
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
},{"./edge.js":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/edge.js","./initializePoints.js":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/initializePoints.js","./utilities.js":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/utilities.js","delaunay-triangulate":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/triangulate.js"}],"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/edge.js":[function(require,module,exports){
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
    }
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
},{"./vector.js":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/vector.js"}],"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/initializeAnts.js":[function(require,module,exports){
'use strict'

var antFunction = require('./ant.js');

// var NBANTS = 4000;

module.exports = function (container, pointsInfos, options) {

	var Ant = antFunction(container, pointsInfos, options);

	var population = new Array(options.nbAnts);
	var possibleStartPointsId = pointsInfos.possibleStartPointsId;

	for (var i = 0; i < options.nbAnts; i++) {
		var newAnt = new Ant(Ant.generateRandStartPoint());
		newAnt.setDirection();
		population[i] = newAnt;
	}

	return population;

}
},{"./ant.js":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/ant.js"}],"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/initializePoints.js":[function(require,module,exports){
'use strict'

var parse = require('parse-svg-path');

var range = require('./utilities.js').range;

var Point = require('./point.js');

var random = Math.random;

// var nbRandomPoints = 500;
// var nbStartPoints = 20;

var nbCity = 2;

var textMesh = true;

// Frame definition
var xInit = 0, yInit = 0;
var w = 1,
    h = 1;

var Achar = "c 4.6011,-11.71047 9.20835,-23.42006 13.8199,-35.12898 4.61156,-11.70892 9.22741,-23.41718 13.84573,-35.125 4.61831,-11.70782 9.23908,-23.41519 13.86046,-35.12233 4.62138,-11.70714 9.24336,-23.41406 13.86411,-35.12097 4.62074,-11.70691 9.24025,-23.41381 13.85667,-35.12092 4.61641,-11.70712 9.22974,-23.41444 13.83814,-35.12218 4.60839,-11.70775 9.21186,-23.41592 13.80854,-35.12474 4.59668,-11.70881 9.18658,-23.418277 13.76785,-35.128603 12.99923,-3.357855 24.30069,-6.821144 33.86893,-4.46411 9.56825,2.357035 17.40328,10.534393 23.46967,30.457833 4.17598,10.62114 8.36031,21.23882 12.54942,31.85453 4.1891,10.61571 8.38298,21.22943 12.57807,31.84266 4.19509,10.61323 8.3914,21.22595 12.58535,31.83965 4.19395,10.61369 8.38555,21.22836 12.57123,31.84549 4.18568,10.61712 8.36544,21.2367 12.53572,31.86021 4.17028,10.6235 8.33108,21.25093 12.47883,31.88377 4.14775,10.63283 8.28245,21.27107 12.40053,31.9162 4.11808,10.64512 8.21955,21.29712 12.30085,31.95749 -5.90896,6.95561 -24.61617,1.11298 -35.59372,3 -16.75248,2.84859 -26.96421,-0.41416 -28.40628,-19 -3.17726,-8.42157 -6.35606,-16.87924 -9.47997,-25.34854 -3.12391,-8.46929 -6.19293,-16.9502 -9.15062,-25.41826 -16.46723,-0.80053 -35.17768,-2.74281 -53.33727,-3.11439 -18.1596,-0.37158 -35.76834,0.82753 -50.03214,6.30976 -4.15679,10.93124 -8.13806,21.9269 -12.15785,32.90546 -4.0198,10.97855 -8.07813,21.94001 -12.38903,32.80282 -9.52758,0.82747 -19.10457,0.96423 -28.69282,0.93363 -9.58824,-0.0306 -19.18773,-0.22855 -28.7603,-0.0705 l 0,-2.19791 z m 174,-109 c -3.17462,-9.51136 -6.44835,-18.99174 -9.7543,-28.46224 -6.81559,-19.51412 -6.24202,-25.37946 -12.3174,-43.60788 -3.16722,-9.51543 -13.60257,-32.32866 -16.49973,-41.92988 -6.08989,12.88576 -11.132,26.59272 -15.93415,40.47476 -4.80215,13.88203 -9.36435,27.93915 -14.49442,41.52524 -2.39425,12.05801 -22.8815,40.11116 2.18745,34 11.05256,-0.44486 22.52939,0.11061 33.85624,0.24955 11.32684,0.13895 22.50369,-0.13863 32.95631,-2.24955 z ";
var Nchar = "c 0,0 0,-23.83334 0,-35.75 0,-11.91667 0,-23.83334 0,-35.75 0,-11.91666 0,-23.83333 0,-35.75 0,-11.916667 0,-23.833335 0,-35.750003 14.64729,1.05313 31.03193,-2.08808 44.60679,1.55415 6.18076,7.610996 12.35438,15.231342 18.51973,22.861013 6.16535,7.62967 12.32244,15.26866 18.47014,22.91696 6.1477,7.64829 12.28602,15.30588 18.41383,22.97274 6.12781,7.66686 12.24512,15.343 18.35081,23.02838 6.10568,7.68538 12.19975,15.38001 18.28107,23.08386 6.08132,7.70384 12.1499,15.41691 18.20462,23.13918 6.05472,7.72226 12.09557,15.45372 18.12145,23.19435 6.02587,7.74063 12.03676,15.49043 18.03156,23.24937 0.4558,-15.4914 0.72655,-30.98725 0.88119,-46.48582 0.15464,-15.49858 0.19319,-30.99988 0.18458,-46.50218 -0.009,-15.5023 -0.0643,-31.0056 -0.0983,-46.50817 -0.034,-15.50258 -0.0461,-31.004432 0.0325,-46.503833 9.33333,0 18.66667,0 28,0 9.33333,0 18.66666,0 28,0 0,11.916667 0,23.833332 0,35.750003 0,11.91666 0,23.83333 0,35.75 0,11.91666 0,23.83333 0,35.75 0,11.91667 0,23.83333 0,35.75 0,11.91666 0,23.83333 0,35.75 0,11.91667 0,23.83334 0,35.75 0,11.91666 0,23.83333 0,35.75 0,11.91667 0,23.83334 0,35.75 -14.6441,-1.05348 -31.026,2.08847 -44.59764,-1.55416 -6.27247,-7.38037 -12.48436,-14.81692 -18.65123,-22.29507 -6.16688,-7.47815 -12.28873,-14.9979 -18.38111,-22.54465 -6.09238,-7.54675 -12.15528,-15.12051 -18.20425,-22.70667 -6.04896,-7.58617 -12.08399,-15.18476 -18.12063,-22.78116 -6.03664,-7.5964 -12.07489,-15.19062 -18.1303,-22.76807 -6.05541,-7.57745 -12.12797,-15.13813 -18.23323,-22.66744 -6.10526,-7.52932 -12.24322,-15.02727 -18.42942,-22.47926 -6.18619,-7.45199 -12.42063,-14.85803 -18.71886,-22.20352 -0.22791,15.16496 -0.36674,30.3308 -0.4489,45.49718 -0.0822,15.16637 -0.10766,30.33329 -0.10889,45.5004 -0.001,15.16712 0.0218,30.33443 0.0367,45.50162 0.0149,15.16718 0.0216,30.33422 -0.0122,45.5008 -9.33333,0 -18.66667,0 -28,0 -9.33333,0 -18.66666,0 -28,0 0,-11.91667 0,-23.83334 0,-35.75 0,-11.91667 0,-23.83334 0,-35.75 0,-11.91666 0,-23.83333 0,-35.75 0,-11.91667 0,-35.75 -10e-6,-35.75 z ";
var Tchar = "c 0,-9.91667 0,-19.83334 0,-29.75 0,-9.91667 0,-19.83334 0,-29.75 0,-9.91666 0,-19.83333 0,-29.75 0,-9.91666 -1.47394,-19.83333 -1.47394,-29.75 0,-15.33334 -29.19273,0 -44.52606,0 -15.33333,0 -30.66667,0 -46,0 0,-8 0,-16 0,-24.000002 0,-8 0,-16.000001 0,-24.000001 9.91666,0 19.83333,0 29.75,0 9.91667,0 19.83333,0 29.75,0 9.91666,0 19.83333,0 29.75,0 9.91667,0 19.83333,0 29.75,0 9.91666,0 19.83333,0 29.75,0 9.91667,0 19.83333,0 29.75,0 9.91666,0 19.83333,0 29.75,0 9.91667,0 19.83333,0 29.75,0 0,8 0,16.000001 0,24.000001 0,8.000002 0,16.000002 0,24.000002 -15,0 -30,0 -45,0 -15,0 -44.88809,-13.52565 -45,1.47394 -0.0735,9.85588 -0.10613,18.2399 -0.11249,28.09923 -0.006,9.85932 0.0135,19.72001 0.045,29.58133 0.0315,9.86132 0.0745,19.72329 0.1143,29.58517 0.0399,9.86188 0.0765,19.72367 0.0954,29.58467 0.0189,9.86099 0.0199,19.72118 -0.0117,29.57984 -0.0315,9.85866 -0.0956,19.7158 -0.20688,29.57069 -0.11131,9.85488 -0.26985,19.70752 -0.49033,29.55719 -0.22047,9.84967 -0.50289,19.69636 -0.86193,29.53937 -7.39624,-0.99964 -18.95658,0.96726 -29.70912,1.82971 -10.75253,0.86245 -24.48556,2.02609 -24.86231,-4.79696 -0.66223,-11.99314 -0.66223,-21.54349 -0.49667,-30.48314 0.16556,-8.93965 0.49667,-17.2686 0.49667,-26.81895 0,-9.55035 0,-19.1007 0,-28.65105 0,-9.55035 0,-19.10069 0,-28.65104 z";
var Schar = "c -16.6587,-2.34387 -33.35995,-6.23515 -49.25137,-12.08831 -15.89142,-5.85316 -30.97305,-13.6682 -44.3926,-23.85957 3.73715,-8.02639 7.96065,-15.79371 12.19101,-23.55228 4.23037,-7.75858 8.46759,-15.50842 12.23219,-23.49984 8.82902,6.69162 18.48235,12.75663 28.67412,17.86344 10.19178,5.10681 20.92205,9.25542 31.90485,12.11421 10.9829,2.8588 22.2183,4.4278 33.4206,4.37539 11.2022,-0.0524 22.3713,-1.72622 33.2212,-5.35304 14.2804,-6.92169 17.3033,-19.64436 13.4946,-31.15028 -3.8087,-11.50591 -14.4489,-21.79507 -27.4946,-23.84972 -10.162,-4.07686 -21.1052,-7.14863 -32.1975,-10.11097 -11.0924,-2.96235 -22.334,-5.81526 -33.0931,-9.45439 -10.7591,-3.63914 -21.03562,-8.06449 -30.19779,-14.17171 -9.16216,-6.10723 -17.20996,-13.89632 -23.51158,-24.26293 -5.1794,-10.97058 -7.48544,-22.77583 -7.30512,-34.52966 0.18031,-11.75382 2.84698,-23.4562 7.61299,-34.22103 4.76601,-10.76483 11.63137,-20.5921 20.20905,-28.595692 8.57769,-8.003592 18.86771,-14.183506 30.48305,-17.653621 12.1821,-4.024145 24.8777,-6.357009 37.6988,-7.121901 12.8211,-0.764892 25.7677,0.03819 38.4517,2.285933 12.6841,2.247744 25.1055,5.940153 36.8765,10.953918 11.7709,5.013764 22.8912,11.348885 32.973,18.882053 -4.4745,7.00368 -8.2581,14.70708 -12.1886,22.23499 -3.9306,7.5279 -8.008,14.88031 -13.0701,21.18198 -15.0728,-11.15158 -33.5009,-20.54491 -52.6114,-24.91726 -19.1105,-4.37235 -38.9034,-3.72372 -56.7058,5.20861 -10.0321,7.02789 -13.3348,18.84695 -11.1561,29.53596 2.1787,10.68902 9.8387,20.24799 21.732,22.75572 10.1335,4.43307 21.0241,7.69407 32.1101,10.69023 11.0861,2.99616 22.3676,5.72747 33.283,9.10117 10.9155,3.3737 21.4648,7.38978 31.0865,12.95547 9.6217,5.56569 18.3157,12.68099 25.5204,22.25313 6.6301,9.71635 10.5587,20.88023 12.0064,32.41231 1.4476,11.53207 0.4142,23.43234 -2.8797,34.62146 -3.2939,11.18912 -8.8484,21.66709 -16.443,30.35457 -7.5946,8.68749 -17.2293,15.58449 -28.6837,19.61166 -13.0091,5.92928 -27.0197,8.63869 -41.2728,9.63608 -14.253,0.99738 -28.7484,0.28274 -42.7272,-0.63608 z ";

var svgString = "m 0,0 " + Achar +"m 126,-31 " + Nchar + "m 376,24 "+ Tchar +"m 250,120 "+ Schar;

 
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
                X += command[5];
                Y += command[6];
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
        var scaleY = 0.2;
        var deltaX = 0.25;
        var deltaY = 0.3;

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

    if (nbRandomPoints !== 0) {
        for (var i = 0; i < nbStartPoints; i++){
            possibleStartPointsId.push(Math.floor(nbRandomPoints * random()));
        }
    } else {
        for (var i = 0; i < nbStartPoints; i++){
            possibleStartPointsId.push(Math.floor(nbPoints * random()));
        }
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

},{"./point.js":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/point.js","./utilities.js":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/utilities.js","parse-svg-path":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/parse-svg-path/index.js"}],"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/mouse.js":[function(require,module,exports){
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

},{}],"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/point.js":[function(require,module,exports){
'use strict'

function Point(x, y) {
    this.id = undefined;                
    this.x = x;
    this.y = y;
    this.nexts = [];
}

module.exports = Point;
},{}],"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/rendering.js":[function(require,module,exports){
'use strict'

var antFunction = require('./ant.js');

var random = Math.random;

var RANDOMMVT = 0.003;
var ANTSIZE = 0.002;

module.exports = function(container, initVar){
	
	if(!container)
		throw new TypeError('Missing container');

	var edges = initVar.edges;
	var population = initVar.population;
	var pointsInfos = initVar.pointsInfos;
	var nbAnts = population.length;

	var animID;

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
	

	function tick() {
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

		// for (i = 0; i < nbAnts; i++) {
		//     population[i].transit();
		// }

		// pheromon evaporation
		edges.forEach(function(edge){
			if(edge.pheromon > 0){
				edge.pheromon -= 0.0001;
			}
		});

		// for (i = 0; i < edges.length; i++) {
			

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
		// for(var i=0; i<population.length; ++i) {
			
		// }

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


	function setAntCount(opts){

		var previousCount = population.length;

		console.log('nbAnts :', opts.nbAnts);
		console.log('previousCount :', previousCount);


		if (opts.nbAnts > previousCount){
			console.log('pointsInfos :', pointsInfos);
			var Ant = antFunction(container, pointsInfos, opts);

			for (var i = 0; i < opts.nbAnts - previousCount; i++) {
				console.log('test');
				var newAnt = new Ant(Ant.generateRandStartPoint());
				newAnt.setDirection();
				population.push(newAnt);
			}
		}
		else{
			population = population.slice(0, opts.nbAnts);
			nbAnts = population.length;
			console.log('Nb Ants :', population.length);
		}
		console.log('Nouveau compte :', population.length);
	}

	function modifyAnts(opts){
		setAntCount(opts);

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
		setAntCount: setAntCount,
		modifyAnts: modifyAnts
	}
}

},{"./ant.js":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/ant.js"}],"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/utilities.js":[function(require,module,exports){
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
},{}],"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/vector.js":[function(require,module,exports){
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
},{}]},{},["/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/index.js"])
//# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJzb3VyY2VzIjpbIi91c3IvbG9jYWwvbGliL25vZGVfbW9kdWxlcy93YXRjaGlmeS9ub2RlX21vZHVsZXMvYnJvd3NlcmlmeS9ub2RlX21vZHVsZXMvYnJvd3Nlci1wYWNrL19wcmVsdWRlLmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9Qcm9qZXRzX0FudHMvQW50Q29sb255L2luZGV4LmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9Qcm9qZXRzX0FudHMvQW50Q29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvaW5jcmVtZW50YWwtY29udmV4LWh1bGwvaWNoLmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9Qcm9qZXRzX0FudHMvQW50Q29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvaW5jcmVtZW50YWwtY29udmV4LWh1bGwvbm9kZV9tb2R1bGVzL3JvYnVzdC1vcmllbnRhdGlvbi9ub2RlX21vZHVsZXMvcm9idXN0LXNjYWxlL25vZGVfbW9kdWxlcy90d28tc3VtL3R3by1zdW0uanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL1Byb2pldHNfQW50cy9BbnRDb2xvbnkvbm9kZV9tb2R1bGVzL2RlbGF1bmF5LXRyaWFuZ3VsYXRlL25vZGVfbW9kdWxlcy9pbmNyZW1lbnRhbC1jb252ZXgtaHVsbC9ub2RlX21vZHVsZXMvcm9idXN0LW9yaWVudGF0aW9uL25vZGVfbW9kdWxlcy9yb2J1c3Qtc2NhbGUvcm9idXN0LXNjYWxlLmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9Qcm9qZXRzX0FudHMvQW50Q29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvaW5jcmVtZW50YWwtY29udmV4LWh1bGwvbm9kZV9tb2R1bGVzL3JvYnVzdC1vcmllbnRhdGlvbi9ub2RlX21vZHVsZXMvcm9idXN0LXN1YnRyYWN0L3JvYnVzdC1kaWZmLmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9Qcm9qZXRzX0FudHMvQW50Q29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvaW5jcmVtZW50YWwtY29udmV4LWh1bGwvbm9kZV9tb2R1bGVzL3JvYnVzdC1vcmllbnRhdGlvbi9ub2RlX21vZHVsZXMvcm9idXN0LXN1bS9yb2J1c3Qtc3VtLmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9Qcm9qZXRzX0FudHMvQW50Q29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvaW5jcmVtZW50YWwtY29udmV4LWh1bGwvbm9kZV9tb2R1bGVzL3JvYnVzdC1vcmllbnRhdGlvbi9ub2RlX21vZHVsZXMvdHdvLXByb2R1Y3QvdHdvLXByb2R1Y3QuanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL1Byb2pldHNfQW50cy9BbnRDb2xvbnkvbm9kZV9tb2R1bGVzL2RlbGF1bmF5LXRyaWFuZ3VsYXRlL25vZGVfbW9kdWxlcy9pbmNyZW1lbnRhbC1jb252ZXgtaHVsbC9ub2RlX21vZHVsZXMvcm9idXN0LW9yaWVudGF0aW9uL29yaWVudGF0aW9uLmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9Qcm9qZXRzX0FudHMvQW50Q29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvaW5jcmVtZW50YWwtY29udmV4LWh1bGwvbm9kZV9tb2R1bGVzL3NpbXBsaWNpYWwtY29tcGxleC9ub2RlX21vZHVsZXMvYml0LXR3aWRkbGUvdHdpZGRsZS5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vUHJvamV0c19BbnRzL0FudENvbG9ueS9ub2RlX21vZHVsZXMvZGVsYXVuYXktdHJpYW5ndWxhdGUvbm9kZV9tb2R1bGVzL2luY3JlbWVudGFsLWNvbnZleC1odWxsL25vZGVfbW9kdWxlcy9zaW1wbGljaWFsLWNvbXBsZXgvbm9kZV9tb2R1bGVzL3VuaW9uLWZpbmQvaW5kZXguanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL1Byb2pldHNfQW50cy9BbnRDb2xvbnkvbm9kZV9tb2R1bGVzL2RlbGF1bmF5LXRyaWFuZ3VsYXRlL25vZGVfbW9kdWxlcy9pbmNyZW1lbnRhbC1jb252ZXgtaHVsbC9ub2RlX21vZHVsZXMvc2ltcGxpY2lhbC1jb21wbGV4L3RvcG9sb2d5LmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9Qcm9qZXRzX0FudHMvQW50Q29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvdW5pcS91bmlxLmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9Qcm9qZXRzX0FudHMvQW50Q29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS90cmlhbmd1bGF0ZS5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vUHJvamV0c19BbnRzL0FudENvbG9ueS9ub2RlX21vZHVsZXMvcGFyc2Utc3ZnLXBhdGgvaW5kZXguanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL1Byb2pldHNfQW50cy9BbnRDb2xvbnkvc3JjL2FudC5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vUHJvamV0c19BbnRzL0FudENvbG9ueS9zcmMvY3JlYXRlRWRnZXMuanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL1Byb2pldHNfQW50cy9BbnRDb2xvbnkvc3JjL2VkZ2UuanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL1Byb2pldHNfQW50cy9BbnRDb2xvbnkvc3JjL2luaXRpYWxpemVBbnRzLmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9Qcm9qZXRzX0FudHMvQW50Q29sb255L3NyYy9pbml0aWFsaXplUG9pbnRzLmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9Qcm9qZXRzX0FudHMvQW50Q29sb255L3NyYy9tb3VzZS5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vUHJvamV0c19BbnRzL0FudENvbG9ueS9zcmMvcG9pbnQuanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL1Byb2pldHNfQW50cy9BbnRDb2xvbnkvc3JjL3JlbmRlcmluZy5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vUHJvamV0c19BbnRzL0FudENvbG9ueS9zcmMvdXRpbGl0aWVzLmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9Qcm9qZXRzX0FudHMvQW50Q29sb255L3NyYy92ZWN0b3IuanMiXSwibmFtZXMiOltdLCJtYXBwaW5ncyI6IkFBQUE7QUNBQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDdENBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDN2JBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDaEJBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDakRBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUMzSkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQzNKQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDaENBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQzdMQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUM1TUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDekRBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ3RWQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUN6REE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQzlKQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ3ZEQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQzdNQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ3JGQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDN0VBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ3JCQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUMzS0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNuQkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDVEE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDek1BO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDNUJBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EiLCJmaWxlIjoiZ2VuZXJhdGVkLmpzIiwic291cmNlUm9vdCI6IiIsInNvdXJjZXNDb250ZW50IjpbIihmdW5jdGlvbiBlKHQsbixyKXtmdW5jdGlvbiBzKG8sdSl7aWYoIW5bb10pe2lmKCF0W29dKXt2YXIgYT10eXBlb2YgcmVxdWlyZT09XCJmdW5jdGlvblwiJiZyZXF1aXJlO2lmKCF1JiZhKXJldHVybiBhKG8sITApO2lmKGkpcmV0dXJuIGkobywhMCk7dmFyIGY9bmV3IEVycm9yKFwiQ2Fubm90IGZpbmQgbW9kdWxlICdcIitvK1wiJ1wiKTt0aHJvdyBmLmNvZGU9XCJNT0RVTEVfTk9UX0ZPVU5EXCIsZn12YXIgbD1uW29dPXtleHBvcnRzOnt9fTt0W29dWzBdLmNhbGwobC5leHBvcnRzLGZ1bmN0aW9uKGUpe3ZhciBuPXRbb11bMV1bZV07cmV0dXJuIHMobj9uOmUpfSxsLGwuZXhwb3J0cyxlLHQsbixyKX1yZXR1cm4gbltvXS5leHBvcnRzfXZhciBpPXR5cGVvZiByZXF1aXJlPT1cImZ1bmN0aW9uXCImJnJlcXVpcmU7Zm9yKHZhciBvPTA7bzxyLmxlbmd0aDtvKyspcyhyW29dKTtyZXR1cm4gc30pIiwiJ3VzZSBzdHJpY3QnO1xuXG52YXIgaW5pdFJlbmRlcmluZyA9IHJlcXVpcmUoJy4vc3JjL3JlbmRlcmluZy5qcycpO1xudmFyIGluaXRpYWxpemVQb2ludHMgPSByZXF1aXJlKCcuL3NyYy9pbml0aWFsaXplUG9pbnRzLmpzJyk7XG52YXIgY3JlYXRlRWRnZXMgPSByZXF1aXJlKCcuL3NyYy9jcmVhdGVFZGdlcy5qcycpO1xudmFyIGluaXRBbnRzID0gcmVxdWlyZSgnLi9zcmMvaW5pdGlhbGl6ZUFudHMnKTtcblxubW9kdWxlLmV4cG9ydHMgPSBmdW5jdGlvbiBpbml0KGNvbnRhaW5lckVsZW1lbnQsIG9wdGlvbnMpe1xuXG5cdHZhciByZW5kZXIsIHBvaW50c0luZm9zLCBlZGdlcywgcG9wdWxhdGlvbiwgaW5pdFZhcjtcblxuXG5cdGZ1bmN0aW9uIF9pbml0KGNvbnRhaW5lckVsZW1lbnQsIG9wdGlvbnMpe1xuXHRcdHBvaW50c0luZm9zID0gaW5pdGlhbGl6ZVBvaW50cyhvcHRpb25zLm5iU3RhcnQsIG9wdGlvbnMubmJSYW5kKTtcblx0XHRlZGdlcyA9IGNyZWF0ZUVkZ2VzKHBvaW50c0luZm9zLnBvaW50cyk7XG5cdFx0cG9wdWxhdGlvbiA9IGluaXRBbnRzKGNvbnRhaW5lckVsZW1lbnQsIHBvaW50c0luZm9zLCBvcHRpb25zKTtcblx0XHRpbml0VmFyID0ge1xuXHRcdFx0cG9pbnRzSW5mb3M6IHBvaW50c0luZm9zLFxuXHRcdFx0ZWRnZXM6IGVkZ2VzLFxuXHRcdFx0cG9wdWxhdGlvbjogcG9wdWxhdGlvblxuXHRcdH07XG5cdFx0cmVuZGVyID0gaW5pdFJlbmRlcmluZyhjb250YWluZXJFbGVtZW50LCBpbml0VmFyKTtcblx0fVxuXG5cdF9pbml0KGNvbnRhaW5lckVsZW1lbnQsIG9wdGlvbnMpO1xuXG5cdHJldHVybiB7XG5cdFx0dG9nZ2xlUGxheVBhdXNlOiBmdW5jdGlvbigpeyByZW5kZXIudG9nZ2xlUGxheVBhdXNlKCkgfSxcblx0XHRjaGFuZ2VPcHRpb25zOiBmdW5jdGlvbihvcHRzKXtcblx0XHRcdHJlbmRlci5tb2RpZnlBbnRzKG9wdHMpO1xuXHRcdH0sXG5cdFx0cmVzZXQ6IGZ1bmN0aW9uKG9wdHMpe1xuXHRcdFx0cmVuZGVyLnJlc2V0KCk7XG5cblx0XHRcdFx0Ly8gcmVzZXQgZWxlbWVudHNcblx0XHRcdF9pbml0KGNvbnRhaW5lckVsZW1lbnQsIG9wdHMpO1xuXHRcdH1cblx0fTtcbn07IiwiXCJ1c2Ugc3RyaWN0XCJcblxuLy9IaWdoIGxldmVsIGlkZWE6XG4vLyAxLiBVc2UgQ2xhcmtzb24ncyBpbmNyZW1lbnRhbCBjb25zdHJ1Y3Rpb24gdG8gZmluZCBjb252ZXggaHVsbFxuLy8gMi4gUG9pbnQgbG9jYXRpb24gaW4gdHJpYW5ndWxhdGlvbiBieSBqdW1wIGFuZCB3YWxrXG5cbm1vZHVsZS5leHBvcnRzID0gaW5jcmVtZW50YWxDb252ZXhIdWxsXG5cbnZhciBvcmllbnQgPSByZXF1aXJlKFwicm9idXN0LW9yaWVudGF0aW9uXCIpXG52YXIgY29tcGFyZUNlbGwgPSByZXF1aXJlKFwic2ltcGxpY2lhbC1jb21wbGV4XCIpLmNvbXBhcmVDZWxsc1xuXG5mdW5jdGlvbiBjb21wYXJlSW50KGEsIGIpIHtcbiAgcmV0dXJuIGEgLSBiXG59XG5cbmZ1bmN0aW9uIFNpbXBsZXgodmVydGljZXMsIGFkamFjZW50LCBib3VuZGFyeSkge1xuICB0aGlzLnZlcnRpY2VzID0gdmVydGljZXNcbiAgdGhpcy5hZGphY2VudCA9IGFkamFjZW50XG4gIHRoaXMuYm91bmRhcnkgPSBib3VuZGFyeVxuICB0aGlzLmxhc3RWaXNpdGVkID0gLTFcbn1cblxuU2ltcGxleC5wcm90b3R5cGUuZmxpcCA9IGZ1bmN0aW9uKCkge1xuICB2YXIgdCA9IHRoaXMudmVydGljZXNbMF1cbiAgdGhpcy52ZXJ0aWNlc1swXSA9IHRoaXMudmVydGljZXNbMV1cbiAgdGhpcy52ZXJ0aWNlc1sxXSA9IHRcbiAgdmFyIHUgPSB0aGlzLmFkamFjZW50WzBdXG4gIHRoaXMuYWRqYWNlbnRbMF0gPSB0aGlzLmFkamFjZW50WzFdXG4gIHRoaXMuYWRqYWNlbnRbMV0gPSB1XG59XG5cbmZ1bmN0aW9uIEdsdWVGYWNldCh2ZXJ0aWNlcywgY2VsbCwgaW5kZXgpIHtcbiAgdGhpcy52ZXJ0aWNlcyA9IHZlcnRpY2VzXG4gIHRoaXMuY2VsbCA9IGNlbGxcbiAgdGhpcy5pbmRleCA9IGluZGV4XG59XG5cbmZ1bmN0aW9uIGNvbXBhcmVHbHVlKGEsIGIpIHtcbiAgcmV0dXJuIGNvbXBhcmVDZWxsKGEudmVydGljZXMsIGIudmVydGljZXMpXG59XG5cbmZ1bmN0aW9uIGJha2VPcmllbnQoZCkge1xuICB2YXIgY29kZSA9IFtcImZ1bmN0aW9uIG9yaWVudCgpe3ZhciB0dXBsZT10aGlzLnR1cGxlO3JldHVybiB0ZXN0KFwiXVxuICBmb3IodmFyIGk9MDsgaTw9ZDsgKytpKSB7XG4gICAgaWYoaSA+IDApIHtcbiAgICAgIGNvZGUucHVzaChcIixcIilcbiAgICB9XG4gICAgY29kZS5wdXNoKFwidHVwbGVbXCIsIGksIFwiXVwiKVxuICB9XG4gIGNvZGUucHVzaChcIil9cmV0dXJuIG9yaWVudFwiKVxuICB2YXIgcHJvYyA9IG5ldyBGdW5jdGlvbihcInRlc3RcIiwgY29kZS5qb2luKFwiXCIpKVxuICB2YXIgdGVzdCA9IG9yaWVudFtkKzFdXG4gIGlmKCF0ZXN0KSB7XG4gICAgdGVzdCA9IG9yaWVudFxuICB9XG4gIHJldHVybiBwcm9jKHRlc3QpXG59XG5cbnZhciBCQUtFRCA9IFtdXG5cbmZ1bmN0aW9uIFRyaWFuZ3VsYXRpb24oZGltZW5zaW9uLCB2ZXJ0aWNlcywgc2ltcGxpY2VzKSB7XG4gIHRoaXMuZGltZW5zaW9uID0gZGltZW5zaW9uXG4gIHRoaXMudmVydGljZXMgPSB2ZXJ0aWNlc1xuICB0aGlzLnNpbXBsaWNlcyA9IHNpbXBsaWNlc1xuICB0aGlzLmludGVyaW9yID0gc2ltcGxpY2VzLmZpbHRlcihmdW5jdGlvbihjKSB7XG4gICAgcmV0dXJuICFjLmJvdW5kYXJ5XG4gIH0pXG5cbiAgdGhpcy50dXBsZSA9IG5ldyBBcnJheShkaW1lbnNpb24rMSlcbiAgZm9yKHZhciBpPTA7IGk8PWRpbWVuc2lvbjsgKytpKSB7XG4gICAgdGhpcy50dXBsZVtpXSA9IHRoaXMudmVydGljZXNbaV1cbiAgfVxuXG4gIHZhciBvID0gQkFLRURbZGltZW5zaW9uXVxuICBpZighbykge1xuICAgIG8gPSBCQUtFRFtkaW1lbnNpb25dID0gYmFrZU9yaWVudChkaW1lbnNpb24pXG4gIH1cbiAgdGhpcy5vcmllbnQgPSBvXG59XG5cbnZhciBwcm90byA9IFRyaWFuZ3VsYXRpb24ucHJvdG90eXBlXG5cbi8vRGVnZW5lcmF0ZSBzaXR1YXRpb24gd2hlcmUgd2UgYXJlIG9uIGJvdW5kYXJ5LCBidXQgY29wbGFuYXIgdG8gZmFjZVxucHJvdG8uaGFuZGxlQm91bmRhcnlEZWdlbmVyYWN5ID0gZnVuY3Rpb24oY2VsbCwgcG9pbnQpIHtcbiAgdmFyIGQgPSB0aGlzLmRpbWVuc2lvblxuICB2YXIgbiA9IHRoaXMudmVydGljZXMubGVuZ3RoIC0gMVxuICB2YXIgdHVwbGUgPSB0aGlzLnR1cGxlXG4gIHZhciB2ZXJ0cyA9IHRoaXMudmVydGljZXNcblxuICAvL0R1bWIgc29sdXRpb246IEp1c3QgZG8gZGZzIGZyb20gYm91bmRhcnkgY2VsbCB1bnRpbCB3ZSBmaW5kIGFueSBwZWFrLCBvciB0ZXJtaW5hdGVcbiAgdmFyIHRvVmlzaXQgPSBbIGNlbGwgXVxuICBjZWxsLmxhc3RWaXNpdGVkID0gLW5cbiAgd2hpbGUodG9WaXNpdC5sZW5ndGggPiAwKSB7XG4gICAgY2VsbCA9IHRvVmlzaXQucG9wKClcbiAgICB2YXIgY2VsbFZlcnRzID0gY2VsbC52ZXJ0aWNlc1xuICAgIHZhciBjZWxsQWRqID0gY2VsbC5hZGphY2VudFxuICAgIGZvcih2YXIgaT0wOyBpPD1kOyArK2kpIHtcbiAgICAgIHZhciBuZWlnaGJvciA9IGNlbGxBZGpbaV1cbiAgICAgIGlmKCFuZWlnaGJvci5ib3VuZGFyeSB8fCBuZWlnaGJvci5sYXN0VmlzaXRlZCA8PSAtbikge1xuICAgICAgICBjb250aW51ZVxuICAgICAgfVxuICAgICAgdmFyIG52ID0gbmVpZ2hib3IudmVydGljZXNcbiAgICAgIGZvcih2YXIgaj0wOyBqPD1kOyArK2opIHtcbiAgICAgICAgdmFyIHZ2ID0gbnZbal1cbiAgICAgICAgaWYodnYgPCAwKSB7XG4gICAgICAgICAgdHVwbGVbal0gPSBwb2ludFxuICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgIHR1cGxlW2pdID0gdmVydHNbdnZdXG4gICAgICAgIH1cbiAgICAgIH1cbiAgICAgIHZhciBvID0gdGhpcy5vcmllbnQoKVxuICAgICAgaWYobyA+IDApIHtcbiAgICAgICAgcmV0dXJuIG5laWdoYm9yXG4gICAgICB9XG4gICAgICBuZWlnaGJvci5sYXN0VmlzaXRlZCA9IC1uXG4gICAgICBpZihvID09PSAwKSB7XG4gICAgICAgIHRvVmlzaXQucHVzaChuZWlnaGJvcilcbiAgICAgIH1cbiAgICB9XG4gIH1cbiAgcmV0dXJuIG51bGxcbn1cblxucHJvdG8ud2FsayA9IGZ1bmN0aW9uKHBvaW50LCByYW5kb20pIHtcbiAgLy9BbGlhcyBsb2NhbCBwcm9wZXJ0aWVzXG4gIHZhciBuID0gdGhpcy52ZXJ0aWNlcy5sZW5ndGggLSAxXG4gIHZhciBkID0gdGhpcy5kaW1lbnNpb25cbiAgdmFyIHZlcnRzID0gdGhpcy52ZXJ0aWNlc1xuICB2YXIgdHVwbGUgPSB0aGlzLnR1cGxlXG5cbiAgLy9Db21wdXRlIGluaXRpYWwganVtcCBjZWxsXG4gIHZhciBpbml0SW5kZXggPSByYW5kb20gPyAodGhpcy5pbnRlcmlvci5sZW5ndGggKiBNYXRoLnJhbmRvbSgpKXwwIDogKHRoaXMuaW50ZXJpb3IubGVuZ3RoLTEpXG4gIHZhciBjZWxsID0gdGhpcy5pbnRlcmlvclsgaW5pdEluZGV4IF1cblxuICAvL1N0YXJ0IHdhbGtpbmdcbm91dGVyTG9vcDpcbiAgd2hpbGUoIWNlbGwuYm91bmRhcnkpIHtcbiAgICB2YXIgY2VsbFZlcnRzID0gY2VsbC52ZXJ0aWNlc1xuICAgIHZhciBjZWxsQWRqID0gY2VsbC5hZGphY2VudFxuXG4gICAgZm9yKHZhciBpPTA7IGk8PWQ7ICsraSkge1xuICAgICAgdHVwbGVbaV0gPSB2ZXJ0c1tjZWxsVmVydHNbaV1dXG4gICAgfVxuICAgIGNlbGwubGFzdFZpc2l0ZWQgPSBuXG5cbiAgICAvL0ZpbmQgZmFydGhlc3QgYWRqYWNlbnQgY2VsbFxuICAgIGZvcih2YXIgaT0wOyBpPD1kOyArK2kpIHtcbiAgICAgIHZhciBuZWlnaGJvciA9IGNlbGxBZGpbaV1cbiAgICAgIGlmKG5laWdoYm9yLmxhc3RWaXNpdGVkID49IG4pIHtcbiAgICAgICAgY29udGludWVcbiAgICAgIH1cbiAgICAgIHZhciBwcmV2ID0gdHVwbGVbaV1cbiAgICAgIHR1cGxlW2ldID0gcG9pbnRcbiAgICAgIHZhciBvID0gdGhpcy5vcmllbnQoKVxuICAgICAgdHVwbGVbaV0gPSBwcmV2XG4gICAgICBpZihvIDwgMCkge1xuICAgICAgICBjZWxsID0gbmVpZ2hib3JcbiAgICAgICAgY29udGludWUgb3V0ZXJMb29wXG4gICAgICB9IGVsc2Uge1xuICAgICAgICBpZighbmVpZ2hib3IuYm91bmRhcnkpIHtcbiAgICAgICAgICBuZWlnaGJvci5sYXN0VmlzaXRlZCA9IG5cbiAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICBuZWlnaGJvci5sYXN0VmlzaXRlZCA9IC1uXG4gICAgICAgIH1cbiAgICAgIH1cbiAgICB9XG4gICAgcmV0dXJuXG4gIH1cblxuICByZXR1cm4gY2VsbFxufVxuXG5wcm90by5hZGRQZWFrcyA9IGZ1bmN0aW9uKHBvaW50LCBjZWxsKSB7XG4gIHZhciBuID0gdGhpcy52ZXJ0aWNlcy5sZW5ndGggLSAxXG4gIHZhciBkID0gdGhpcy5kaW1lbnNpb25cbiAgdmFyIHZlcnRzID0gdGhpcy52ZXJ0aWNlc1xuICB2YXIgdHVwbGUgPSB0aGlzLnR1cGxlXG4gIHZhciBpbnRlcmlvciA9IHRoaXMuaW50ZXJpb3JcbiAgdmFyIHNpbXBsaWNlcyA9IHRoaXMuc2ltcGxpY2VzXG5cbiAgLy9XYWxraW5nIGZpbmlzaGVkIGF0IGJvdW5kYXJ5LCB0aW1lIHRvIGFkZCBwZWFrc1xuICB2YXIgdG92aXNpdCA9IFsgY2VsbCBdXG5cbiAgLy9TdHJldGNoIGluaXRpYWwgYm91bmRhcnkgY2VsbCBpbnRvIGEgcGVha1xuICBjZWxsLmxhc3RWaXNpdGVkID0gblxuICBjZWxsLnZlcnRpY2VzW2NlbGwudmVydGljZXMuaW5kZXhPZigtMSldID0gblxuICBjZWxsLmJvdW5kYXJ5ID0gZmFsc2VcbiAgaW50ZXJpb3IucHVzaChjZWxsKVxuXG4gIC8vUmVjb3JkIGEgbGlzdCBvZiBhbGwgbmV3IGJvdW5kYXJpZXMgY3JlYXRlZCBieSBhZGRlZCBwZWFrcyBzbyB3ZSBjYW4gZ2x1ZSB0aGVtIHRvZ2V0aGVyIHdoZW4gd2UgYXJlIGFsbCBkb25lXG4gIHZhciBnbHVlRmFjZXRzID0gW11cblxuICAvL0RvIGEgdHJhdmVyc2FsIG9mIHRoZSBib3VuZGFyeSB3YWxraW5nIG91dHdhcmQgZnJvbSBzdGFydGluZyBwZWFrXG4gIHdoaWxlKHRvdmlzaXQubGVuZ3RoID4gMCkge1xuICAgIC8vUG9wIG9mZiBwZWFrIGFuZCB3YWxrIG92ZXIgYWRqYWNlbnQgY2VsbHNcbiAgICB2YXIgY2VsbCA9IHRvdmlzaXQucG9wKClcbiAgICB2YXIgY2VsbFZlcnRzID0gY2VsbC52ZXJ0aWNlc1xuICAgIHZhciBjZWxsQWRqID0gY2VsbC5hZGphY2VudFxuICAgIHZhciBpbmRleE9mTiA9IGNlbGxWZXJ0cy5pbmRleE9mKG4pXG4gICAgaWYoaW5kZXhPZk4gPCAwKSB7XG4gICAgICBjb250aW51ZVxuICAgIH1cblxuICAgIGZvcih2YXIgaT0wOyBpPD1kOyArK2kpIHtcbiAgICAgIGlmKGkgPT09IGluZGV4T2ZOKSB7XG4gICAgICAgIGNvbnRpbnVlXG4gICAgICB9XG5cbiAgICAgIC8vRm9yIGVhY2ggYm91bmRhcnkgbmVpZ2hib3Igb2YgdGhlIGNlbGxcbiAgICAgIHZhciBuZWlnaGJvciA9IGNlbGxBZGpbaV1cbiAgICAgIGlmKCFuZWlnaGJvci5ib3VuZGFyeSB8fCBuZWlnaGJvci5sYXN0VmlzaXRlZCA+PSBuKSB7XG4gICAgICAgIGNvbnRpbnVlXG4gICAgICB9XG5cbiAgICAgIHZhciBudiA9IG5laWdoYm9yLnZlcnRpY2VzXG5cbiAgICAgIC8vVGVzdCBpZiBuZWlnaGJvciBpcyBhIHBlYWtcbiAgICAgIGlmKG5laWdoYm9yLmxhc3RWaXNpdGVkICE9PSAtbikgeyAgICAgIFxuICAgICAgICAvL0NvbXB1dGUgb3JpZW50YXRpb24gb2YgcCByZWxhdGl2ZSB0byBlYWNoIGJvdW5kYXJ5IHBlYWtcbiAgICAgICAgdmFyIGluZGV4T2ZOZWcxID0gMFxuICAgICAgICBmb3IodmFyIGo9MDsgajw9ZDsgKytqKSB7XG4gICAgICAgICAgaWYobnZbal0gPCAwKSB7XG4gICAgICAgICAgICBpbmRleE9mTmVnMSA9IGpcbiAgICAgICAgICAgIHR1cGxlW2pdID0gcG9pbnRcbiAgICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgdHVwbGVbal0gPSB2ZXJ0c1tudltqXV1cbiAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICAgICAgdmFyIG8gPSB0aGlzLm9yaWVudCgpXG5cbiAgICAgICAgLy9UZXN0IGlmIG5laWdoYm9yIGNlbGwgaXMgYWxzbyBhIHBlYWtcbiAgICAgICAgaWYobyA+IDApIHtcbiAgICAgICAgICBudltpbmRleE9mTmVnMV0gPSBuXG4gICAgICAgICAgbmVpZ2hib3IuYm91bmRhcnkgPSBmYWxzZVxuICAgICAgICAgIGludGVyaW9yLnB1c2gobmVpZ2hib3IpXG4gICAgICAgICAgdG92aXNpdC5wdXNoKG5laWdoYm9yKVxuICAgICAgICAgIG5laWdoYm9yLmxhc3RWaXNpdGVkID0gblxuICAgICAgICAgIGNvbnRpbnVlXG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgbmVpZ2hib3IubGFzdFZpc2l0ZWQgPSAtblxuICAgICAgICB9XG4gICAgICB9XG5cbiAgICAgIHZhciBuYSA9IG5laWdoYm9yLmFkamFjZW50XG5cbiAgICAgIC8vT3RoZXJ3aXNlLCByZXBsYWNlIG5laWdoYm9yIHdpdGggbmV3IGZhY2VcbiAgICAgIHZhciB2dmVydHMgPSBjZWxsVmVydHMuc2xpY2UoKVxuICAgICAgdmFyIHZhZGogPSBjZWxsQWRqLnNsaWNlKClcbiAgICAgIHZhciBuY2VsbCA9IG5ldyBTaW1wbGV4KHZ2ZXJ0cywgdmFkaiwgdHJ1ZSlcbiAgICAgIHNpbXBsaWNlcy5wdXNoKG5jZWxsKVxuXG4gICAgICAvL0Nvbm5lY3QgdG8gbmVpZ2hib3JcbiAgICAgIHZhciBvcHBvc2l0ZSA9IG5hLmluZGV4T2YoY2VsbClcbiAgICAgIGlmKG9wcG9zaXRlIDwgMCkge1xuICAgICAgICBjb250aW51ZVxuICAgICAgfVxuICAgICAgbmFbb3Bwb3NpdGVdID0gbmNlbGxcbiAgICAgIHZhZGpbaW5kZXhPZk5dID0gbmVpZ2hib3JcblxuICAgICAgLy9Db25uZWN0IHRvIGNlbGxcbiAgICAgIHZ2ZXJ0c1tpXSA9IC0xXG4gICAgICB2YWRqW2ldID0gY2VsbFxuICAgICAgY2VsbEFkaltpXSA9IG5jZWxsXG5cbiAgICAgIC8vRmxpcCBmYWNldFxuICAgICAgbmNlbGwuZmxpcCgpXG5cbiAgICAgIC8vQWRkIHRvIGdsdWUgbGlzdFxuICAgICAgZm9yKHZhciBqPTA7IGo8PWQ7ICsraikge1xuICAgICAgICB2YXIgdXUgPSB2dmVydHNbal1cbiAgICAgICAgaWYodXUgPCAwIHx8IHV1ID09PSBuKSB7XG4gICAgICAgICAgY29udGludWVcbiAgICAgICAgfVxuICAgICAgICB2YXIgbmZhY2UgPSBuZXcgQXJyYXkoZC0xKVxuICAgICAgICB2YXIgbnB0ciA9IDBcbiAgICAgICAgZm9yKHZhciBrPTA7IGs8PWQ7ICsraykge1xuICAgICAgICAgIHZhciB2diA9IHZ2ZXJ0c1trXVxuICAgICAgICAgIGlmKHZ2IDwgMCB8fCBrID09PSBqKSB7XG4gICAgICAgICAgICBjb250aW51ZVxuICAgICAgICAgIH1cbiAgICAgICAgICBuZmFjZVtucHRyKytdID0gdnZcbiAgICAgICAgfVxuICAgICAgICBnbHVlRmFjZXRzLnB1c2gobmV3IEdsdWVGYWNldChuZmFjZSwgbmNlbGwsIGopKVxuICAgICAgfVxuICAgIH1cbiAgfVxuXG4gIC8vR2x1ZSBib3VuZGFyeSBmYWNldHMgdG9nZXRoZXJcbiAgZ2x1ZUZhY2V0cy5zb3J0KGNvbXBhcmVHbHVlKVxuXG4gIGZvcih2YXIgaT0wOyBpKzE8Z2x1ZUZhY2V0cy5sZW5ndGg7IGkrPTIpIHtcbiAgICB2YXIgYSA9IGdsdWVGYWNldHNbaV1cbiAgICB2YXIgYiA9IGdsdWVGYWNldHNbaSsxXVxuICAgIHZhciBhaSA9IGEuaW5kZXhcbiAgICB2YXIgYmkgPSBiLmluZGV4XG4gICAgaWYoYWkgPCAwIHx8IGJpIDwgMCkge1xuICAgICAgY29udGludWVcbiAgICB9XG4gICAgYS5jZWxsLmFkamFjZW50W2EuaW5kZXhdID0gYi5jZWxsXG4gICAgYi5jZWxsLmFkamFjZW50W2IuaW5kZXhdID0gYS5jZWxsXG4gIH1cbn1cblxucHJvdG8uaW5zZXJ0ID0gZnVuY3Rpb24ocG9pbnQsIHJhbmRvbSkge1xuICAvL0FkZCBwb2ludFxuICB2YXIgdmVydHMgPSB0aGlzLnZlcnRpY2VzXG4gIHZlcnRzLnB1c2gocG9pbnQpXG5cbiAgdmFyIGNlbGwgPSB0aGlzLndhbGsocG9pbnQsIHJhbmRvbSlcbiAgaWYoIWNlbGwpIHtcbiAgICByZXR1cm5cbiAgfVxuXG4gIC8vQWxpYXMgbG9jYWwgcHJvcGVydGllc1xuICB2YXIgZCA9IHRoaXMuZGltZW5zaW9uXG4gIHZhciB0dXBsZSA9IHRoaXMudHVwbGVcblxuICAvL0RlZ2VuZXJhdGUgY2FzZTogSWYgcG9pbnQgaXMgY29wbGFuYXIgdG8gY2VsbCwgdGhlbiB3YWxrIHVudGlsIHdlIGZpbmQgYSBub24tZGVnZW5lcmF0ZSBib3VuZGFyeVxuICBmb3IodmFyIGk9MDsgaTw9ZDsgKytpKSB7XG4gICAgdmFyIHZ2ID0gY2VsbC52ZXJ0aWNlc1tpXVxuICAgIGlmKHZ2IDwgMCkge1xuICAgICAgdHVwbGVbaV0gPSBwb2ludFxuICAgIH0gZWxzZSB7XG4gICAgICB0dXBsZVtpXSA9IHZlcnRzW3Z2XVxuICAgIH1cbiAgfVxuICB2YXIgbyA9IHRoaXMub3JpZW50KHR1cGxlKVxuICBpZihvIDwgMCkge1xuICAgIHJldHVyblxuICB9IGVsc2UgaWYobyA9PT0gMCkge1xuICAgIGNlbGwgPSB0aGlzLmhhbmRsZUJvdW5kYXJ5RGVnZW5lcmFjeShjZWxsLCBwb2ludClcbiAgICBpZighY2VsbCkge1xuICAgICAgcmV0dXJuXG4gICAgfVxuICB9XG5cbiAgLy9BZGQgcGVha3NcbiAgdGhpcy5hZGRQZWFrcyhwb2ludCwgY2VsbClcbn1cblxuLy9FeHRyYWN0IGFsbCBib3VuZGFyeSBjZWxsc1xucHJvdG8uYm91bmRhcnkgPSBmdW5jdGlvbigpIHtcbiAgdmFyIGQgPSB0aGlzLmRpbWVuc2lvblxuICB2YXIgYm91bmRhcnkgPSBbXVxuICB2YXIgY2VsbHMgPSB0aGlzLnNpbXBsaWNlc1xuICB2YXIgbmMgPSBjZWxscy5sZW5ndGhcbiAgZm9yKHZhciBpPTA7IGk8bmM7ICsraSkge1xuICAgIHZhciBjID0gY2VsbHNbaV1cbiAgICBpZihjLmJvdW5kYXJ5KSB7XG4gICAgICB2YXIgYmNlbGwgPSBuZXcgQXJyYXkoZClcbiAgICAgIHZhciBjdiA9IGMudmVydGljZXNcbiAgICAgIHZhciBwdHIgPSAwXG4gICAgICB2YXIgcGFyaXR5ID0gMFxuICAgICAgZm9yKHZhciBqPTA7IGo8PWQ7ICsraikge1xuICAgICAgICBpZihjdltqXSA+PSAwKSB7XG4gICAgICAgICAgYmNlbGxbcHRyKytdID0gY3Zbal1cbiAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICBwYXJpdHkgPSBqJjFcbiAgICAgICAgfVxuICAgICAgfVxuICAgICAgaWYocGFyaXR5ID09PSAoZCYxKSkge1xuICAgICAgICB2YXIgdCA9IGJjZWxsWzBdXG4gICAgICAgIGJjZWxsWzBdID0gYmNlbGxbMV1cbiAgICAgICAgYmNlbGxbMV0gPSB0XG4gICAgICB9XG4gICAgICBib3VuZGFyeS5wdXNoKGJjZWxsKVxuICAgIH1cbiAgfVxuICByZXR1cm4gYm91bmRhcnlcbn1cblxuZnVuY3Rpb24gaW5jcmVtZW50YWxDb252ZXhIdWxsKHBvaW50cywgcmFuZG9tU2VhcmNoKSB7XG4gIHZhciBuID0gcG9pbnRzLmxlbmd0aFxuICBpZihuID09PSAwKSB7XG4gICAgdGhyb3cgbmV3IEVycm9yKFwiTXVzdCBoYXZlIGF0IGxlYXN0IGQrMSBwb2ludHNcIilcbiAgfVxuICB2YXIgZCA9IHBvaW50c1swXS5sZW5ndGhcbiAgaWYobiA8PSBkKSB7XG4gICAgdGhyb3cgbmV3IEVycm9yKFwiTXVzdCBpbnB1dCBhdCBsZWFzdCBkKzEgcG9pbnRzXCIpXG4gIH1cblxuICAvL0ZJWE1FOiBUaGlzIGNvdWxkIGJlIGRlZ2VuZXJhdGUsIGJ1dCBuZWVkIHRvIHNlbGVjdCBkKzEgbm9uLWNvcGxhbmFyIHBvaW50cyB0byBib290c3RyYXAgcHJvY2Vzc1xuICB2YXIgaW5pdGlhbFNpbXBsZXggPSBwb2ludHMuc2xpY2UoMCwgZCsxKVxuXG4gIC8vTWFrZSBzdXJlIGluaXRpYWwgc2ltcGxleCBpcyBwb3NpdGl2ZWx5IG9yaWVudGVkXG4gIHZhciBvID0gb3JpZW50LmFwcGx5KHZvaWQgMCwgaW5pdGlhbFNpbXBsZXgpXG4gIGlmKG8gPT09IDApIHtcbiAgICB0aHJvdyBuZXcgRXJyb3IoXCJJbnB1dCBub3QgaW4gZ2VuZXJhbCBwb3NpdGlvblwiKVxuICB9XG4gIHZhciBpbml0aWFsQ29vcmRzID0gbmV3IEFycmF5KGQrMSlcbiAgZm9yKHZhciBpPTA7IGk8PWQ7ICsraSkge1xuICAgIGluaXRpYWxDb29yZHNbaV0gPSBpXG4gIH1cbiAgaWYobyA8IDApIHtcbiAgICBpbml0aWFsQ29vcmRzWzBdID0gMVxuICAgIGluaXRpYWxDb29yZHNbMV0gPSAwXG4gIH1cblxuICAvL0NyZWF0ZSBpbml0aWFsIHRvcG9sb2dpY2FsIGluZGV4LCBnbHVlIHBvaW50ZXJzIHRvZ2V0aGVyIChraW5kIG9mIG1lc3N5KVxuICB2YXIgaW5pdGlhbENlbGwgPSBuZXcgU2ltcGxleChpbml0aWFsQ29vcmRzLCBuZXcgQXJyYXkoZCsxKSwgZmFsc2UpXG4gIHZhciBib3VuZGFyeSA9IGluaXRpYWxDZWxsLmFkamFjZW50XG4gIHZhciBsaXN0ID0gbmV3IEFycmF5KGQrMilcbiAgZm9yKHZhciBpPTA7IGk8PWQ7ICsraSkge1xuICAgIHZhciB2ZXJ0cyA9IGluaXRpYWxDb29yZHMuc2xpY2UoKVxuICAgIGZvcih2YXIgaj0wOyBqPD1kOyArK2opIHtcbiAgICAgIGlmKGogPT09IGkpIHtcbiAgICAgICAgdmVydHNbal0gPSAtMVxuICAgICAgfVxuICAgIH1cbiAgICB2YXIgdCA9IHZlcnRzWzBdXG4gICAgdmVydHNbMF0gPSB2ZXJ0c1sxXVxuICAgIHZlcnRzWzFdID0gdFxuICAgIHZhciBjZWxsID0gbmV3IFNpbXBsZXgodmVydHMsIG5ldyBBcnJheShkKzEpLCB0cnVlKVxuICAgIGJvdW5kYXJ5W2ldID0gY2VsbFxuICAgIGxpc3RbaV0gPSBjZWxsXG4gIH1cbiAgbGlzdFtkKzFdID0gaW5pdGlhbENlbGxcbiAgZm9yKHZhciBpPTA7IGk8PWQ7ICsraSkge1xuICAgIHZhciB2ZXJ0cyA9IGJvdW5kYXJ5W2ldLnZlcnRpY2VzXG4gICAgdmFyIGFkaiA9IGJvdW5kYXJ5W2ldLmFkamFjZW50XG4gICAgZm9yKHZhciBqPTA7IGo8PWQ7ICsraikge1xuICAgICAgdmFyIHYgPSB2ZXJ0c1tqXVxuICAgICAgaWYodiA8IDApIHtcbiAgICAgICAgYWRqW2pdID0gaW5pdGlhbENlbGxcbiAgICAgICAgY29udGludWVcbiAgICAgIH1cbiAgICAgIGZvcih2YXIgaz0wOyBrPD1kOyArK2spIHtcbiAgICAgICAgaWYoYm91bmRhcnlba10udmVydGljZXMuaW5kZXhPZih2KSA8IDApIHtcbiAgICAgICAgICBhZGpbal0gPSBib3VuZGFyeVtrXVxuICAgICAgICB9XG4gICAgICB9XG4gICAgfVxuICB9XG5cbiAgLy9Jbml0aWFsaXplIHRyaWFuZ2xlc1xuICB2YXIgdHJpYW5nbGVzID0gbmV3IFRyaWFuZ3VsYXRpb24oZCwgaW5pdGlhbFNpbXBsZXgsIGxpc3QpXG5cbiAgLy9JbnNlcnQgcmVtYWluaW5nIHBvaW50c1xuICB2YXIgdXNlUmFuZG9tID0gISFyYW5kb21TZWFyY2hcbiAgZm9yKHZhciBpPWQrMTsgaTxuOyArK2kpIHtcbiAgICB0cmlhbmdsZXMuaW5zZXJ0KHBvaW50c1tpXSwgdXNlUmFuZG9tKVxuICB9XG4gIFxuICAvL0V4dHJhY3QgYm91bmRhcnkgY2VsbHNcbiAgcmV0dXJuIHRyaWFuZ2xlcy5ib3VuZGFyeSgpXG59IiwiXCJ1c2Ugc3RyaWN0XCJcblxubW9kdWxlLmV4cG9ydHMgPSBmYXN0VHdvU3VtXG5cbmZ1bmN0aW9uIGZhc3RUd29TdW0oYSwgYiwgcmVzdWx0KSB7XG5cdHZhciB4ID0gYSArIGJcblx0dmFyIGJ2ID0geCAtIGFcblx0dmFyIGF2ID0geCAtIGJ2XG5cdHZhciBiciA9IGIgLSBidlxuXHR2YXIgYXIgPSBhIC0gYXZcblx0aWYocmVzdWx0KSB7XG5cdFx0cmVzdWx0WzBdID0gYXIgKyBiclxuXHRcdHJlc3VsdFsxXSA9IHhcblx0XHRyZXR1cm4gcmVzdWx0XG5cdH1cblx0cmV0dXJuIFthciticiwgeF1cbn0iLCJcInVzZSBzdHJpY3RcIlxuXG52YXIgdHdvUHJvZHVjdCA9IHJlcXVpcmUoXCJ0d28tcHJvZHVjdFwiKVxudmFyIHR3b1N1bSA9IHJlcXVpcmUoXCJ0d28tc3VtXCIpXG5cbm1vZHVsZS5leHBvcnRzID0gc2NhbGVMaW5lYXJFeHBhbnNpb25cblxuZnVuY3Rpb24gc2NhbGVMaW5lYXJFeHBhbnNpb24oZSwgc2NhbGUpIHtcbiAgdmFyIG4gPSBlLmxlbmd0aFxuICBpZihuID09PSAxKSB7XG4gICAgdmFyIHRzID0gdHdvUHJvZHVjdChlWzBdLCBzY2FsZSlcbiAgICBpZih0c1swXSkge1xuICAgICAgcmV0dXJuIHRzXG4gICAgfVxuICAgIHJldHVybiBbIHRzWzFdIF1cbiAgfVxuICB2YXIgZyA9IG5ldyBBcnJheSgyICogbilcbiAgdmFyIHEgPSBbMC4xLCAwLjFdXG4gIHZhciB0ID0gWzAuMSwgMC4xXVxuICB2YXIgY291bnQgPSAwXG4gIHR3b1Byb2R1Y3QoZVswXSwgc2NhbGUsIHEpXG4gIGlmKHFbMF0pIHtcbiAgICBnW2NvdW50KytdID0gcVswXVxuICB9XG4gIGZvcih2YXIgaT0xOyBpPG47ICsraSkge1xuICAgIHR3b1Byb2R1Y3QoZVtpXSwgc2NhbGUsIHQpXG4gICAgdmFyIHBxID0gcVsxXVxuICAgIHR3b1N1bShwcSwgdFswXSwgcSlcbiAgICBpZihxWzBdKSB7XG4gICAgICBnW2NvdW50KytdID0gcVswXVxuICAgIH1cbiAgICB2YXIgYSA9IHRbMV1cbiAgICB2YXIgYiA9IHFbMV1cbiAgICB2YXIgeCA9IGEgKyBiXG4gICAgdmFyIGJ2ID0geCAtIGFcbiAgICB2YXIgeSA9IGIgLSBidlxuICAgIHFbMV0gPSB4XG4gICAgaWYoeSkge1xuICAgICAgZ1tjb3VudCsrXSA9IHlcbiAgICB9XG4gIH1cbiAgaWYocVsxXSkge1xuICAgIGdbY291bnQrK10gPSBxWzFdXG4gIH1cbiAgaWYoY291bnQgPT09IDApIHtcbiAgICBnW2NvdW50KytdID0gMC4wXG4gIH1cbiAgZy5sZW5ndGggPSBjb3VudFxuICByZXR1cm4gZ1xufSIsIlwidXNlIHN0cmljdFwiXG5cbm1vZHVsZS5leHBvcnRzID0gcm9idXN0U3VidHJhY3RcblxuLy9FYXN5IGNhc2U6IEFkZCB0d28gc2NhbGFyc1xuZnVuY3Rpb24gc2NhbGFyU2NhbGFyKGEsIGIpIHtcbiAgdmFyIHggPSBhICsgYlxuICB2YXIgYnYgPSB4IC0gYVxuICB2YXIgYXYgPSB4IC0gYnZcbiAgdmFyIGJyID0gYiAtIGJ2XG4gIHZhciBhciA9IGEgLSBhdlxuICB2YXIgeSA9IGFyICsgYnJcbiAgaWYoeSkge1xuICAgIHJldHVybiBbeSwgeF1cbiAgfVxuICByZXR1cm4gW3hdXG59XG5cbmZ1bmN0aW9uIHJvYnVzdFN1YnRyYWN0KGUsIGYpIHtcbiAgdmFyIG5lID0gZS5sZW5ndGh8MFxuICB2YXIgbmYgPSBmLmxlbmd0aHwwXG4gIGlmKG5lID09PSAxICYmIG5mID09PSAxKSB7XG4gICAgcmV0dXJuIHNjYWxhclNjYWxhcihlWzBdLCAtZlswXSlcbiAgfVxuICB2YXIgbiA9IG5lICsgbmZcbiAgdmFyIGcgPSBuZXcgQXJyYXkobilcbiAgdmFyIGNvdW50ID0gMFxuICB2YXIgZXB0ciA9IDBcbiAgdmFyIGZwdHIgPSAwXG4gIHZhciBhYnMgPSBNYXRoLmFic1xuICB2YXIgZWkgPSBlW2VwdHJdXG4gIHZhciBlYSA9IGFicyhlaSlcbiAgdmFyIGZpID0gLWZbZnB0cl1cbiAgdmFyIGZhID0gYWJzKGZpKVxuICB2YXIgYSwgYlxuICBpZihlYSA8IGZhKSB7XG4gICAgYiA9IGVpXG4gICAgZXB0ciArPSAxXG4gICAgaWYoZXB0ciA8IG5lKSB7XG4gICAgICBlaSA9IGVbZXB0cl1cbiAgICAgIGVhID0gYWJzKGVpKVxuICAgIH1cbiAgfSBlbHNlIHtcbiAgICBiID0gZmlcbiAgICBmcHRyICs9IDFcbiAgICBpZihmcHRyIDwgbmYpIHtcbiAgICAgIGZpID0gLWZbZnB0cl1cbiAgICAgIGZhID0gYWJzKGZpKVxuICAgIH1cbiAgfVxuICBpZigoZXB0ciA8IG5lICYmIGVhIDwgZmEpIHx8IChmcHRyID49IG5mKSkge1xuICAgIGEgPSBlaVxuICAgIGVwdHIgKz0gMVxuICAgIGlmKGVwdHIgPCBuZSkge1xuICAgICAgZWkgPSBlW2VwdHJdXG4gICAgICBlYSA9IGFicyhlaSlcbiAgICB9XG4gIH0gZWxzZSB7XG4gICAgYSA9IGZpXG4gICAgZnB0ciArPSAxXG4gICAgaWYoZnB0ciA8IG5mKSB7XG4gICAgICBmaSA9IC1mW2ZwdHJdXG4gICAgICBmYSA9IGFicyhmaSlcbiAgICB9XG4gIH1cbiAgdmFyIHggPSBhICsgYlxuICB2YXIgYnYgPSB4IC0gYVxuICB2YXIgeSA9IGIgLSBidlxuICB2YXIgcTAgPSB5XG4gIHZhciBxMSA9IHhcbiAgdmFyIF94LCBfYnYsIF9hdiwgX2JyLCBfYXJcbiAgd2hpbGUoZXB0ciA8IG5lICYmIGZwdHIgPCBuZikge1xuICAgIGlmKGVhIDwgZmEpIHtcbiAgICAgIGEgPSBlaVxuICAgICAgZXB0ciArPSAxXG4gICAgICBpZihlcHRyIDwgbmUpIHtcbiAgICAgICAgZWkgPSBlW2VwdHJdXG4gICAgICAgIGVhID0gYWJzKGVpKVxuICAgICAgfVxuICAgIH0gZWxzZSB7XG4gICAgICBhID0gZmlcbiAgICAgIGZwdHIgKz0gMVxuICAgICAgaWYoZnB0ciA8IG5mKSB7XG4gICAgICAgIGZpID0gLWZbZnB0cl1cbiAgICAgICAgZmEgPSBhYnMoZmkpXG4gICAgICB9XG4gICAgfVxuICAgIGIgPSBxMFxuICAgIHggPSBhICsgYlxuICAgIGJ2ID0geCAtIGFcbiAgICB5ID0gYiAtIGJ2XG4gICAgaWYoeSkge1xuICAgICAgZ1tjb3VudCsrXSA9IHlcbiAgICB9XG4gICAgX3ggPSBxMSArIHhcbiAgICBfYnYgPSBfeCAtIHExXG4gICAgX2F2ID0gX3ggLSBfYnZcbiAgICBfYnIgPSB4IC0gX2J2XG4gICAgX2FyID0gcTEgLSBfYXZcbiAgICBxMCA9IF9hciArIF9iclxuICAgIHExID0gX3hcbiAgfVxuICB3aGlsZShlcHRyIDwgbmUpIHtcbiAgICBhID0gZWlcbiAgICBiID0gcTBcbiAgICB4ID0gYSArIGJcbiAgICBidiA9IHggLSBhXG4gICAgeSA9IGIgLSBidlxuICAgIGlmKHkpIHtcbiAgICAgIGdbY291bnQrK10gPSB5XG4gICAgfVxuICAgIF94ID0gcTEgKyB4XG4gICAgX2J2ID0gX3ggLSBxMVxuICAgIF9hdiA9IF94IC0gX2J2XG4gICAgX2JyID0geCAtIF9idlxuICAgIF9hciA9IHExIC0gX2F2XG4gICAgcTAgPSBfYXIgKyBfYnJcbiAgICBxMSA9IF94XG4gICAgZXB0ciArPSAxXG4gICAgaWYoZXB0ciA8IG5lKSB7XG4gICAgICBlaSA9IGVbZXB0cl1cbiAgICB9XG4gIH1cbiAgd2hpbGUoZnB0ciA8IG5mKSB7XG4gICAgYSA9IGZpXG4gICAgYiA9IHEwXG4gICAgeCA9IGEgKyBiXG4gICAgYnYgPSB4IC0gYVxuICAgIHkgPSBiIC0gYnZcbiAgICBpZih5KSB7XG4gICAgICBnW2NvdW50KytdID0geVxuICAgIH0gXG4gICAgX3ggPSBxMSArIHhcbiAgICBfYnYgPSBfeCAtIHExXG4gICAgX2F2ID0gX3ggLSBfYnZcbiAgICBfYnIgPSB4IC0gX2J2XG4gICAgX2FyID0gcTEgLSBfYXZcbiAgICBxMCA9IF9hciArIF9iclxuICAgIHExID0gX3hcbiAgICBmcHRyICs9IDFcbiAgICBpZihmcHRyIDwgbmYpIHtcbiAgICAgIGZpID0gLWZbZnB0cl1cbiAgICB9XG4gIH1cbiAgaWYocTApIHtcbiAgICBnW2NvdW50KytdID0gcTBcbiAgfVxuICBpZihxMSkge1xuICAgIGdbY291bnQrK10gPSBxMVxuICB9XG4gIGlmKCFjb3VudCkge1xuICAgIGdbY291bnQrK10gPSAwLjAgIFxuICB9XG4gIGcubGVuZ3RoID0gY291bnRcbiAgcmV0dXJuIGdcbn0iLCJcInVzZSBzdHJpY3RcIlxuXG5tb2R1bGUuZXhwb3J0cyA9IGxpbmVhckV4cGFuc2lvblN1bVxuXG4vL0Vhc3kgY2FzZTogQWRkIHR3byBzY2FsYXJzXG5mdW5jdGlvbiBzY2FsYXJTY2FsYXIoYSwgYikge1xuICB2YXIgeCA9IGEgKyBiXG4gIHZhciBidiA9IHggLSBhXG4gIHZhciBhdiA9IHggLSBidlxuICB2YXIgYnIgPSBiIC0gYnZcbiAgdmFyIGFyID0gYSAtIGF2XG4gIHZhciB5ID0gYXIgKyBiclxuICBpZih5KSB7XG4gICAgcmV0dXJuIFt5LCB4XVxuICB9XG4gIHJldHVybiBbeF1cbn1cblxuZnVuY3Rpb24gbGluZWFyRXhwYW5zaW9uU3VtKGUsIGYpIHtcbiAgdmFyIG5lID0gZS5sZW5ndGh8MFxuICB2YXIgbmYgPSBmLmxlbmd0aHwwXG4gIGlmKG5lID09PSAxICYmIG5mID09PSAxKSB7XG4gICAgcmV0dXJuIHNjYWxhclNjYWxhcihlWzBdLCBmWzBdKVxuICB9XG4gIHZhciBuID0gbmUgKyBuZlxuICB2YXIgZyA9IG5ldyBBcnJheShuKVxuICB2YXIgY291bnQgPSAwXG4gIHZhciBlcHRyID0gMFxuICB2YXIgZnB0ciA9IDBcbiAgdmFyIGFicyA9IE1hdGguYWJzXG4gIHZhciBlaSA9IGVbZXB0cl1cbiAgdmFyIGVhID0gYWJzKGVpKVxuICB2YXIgZmkgPSBmW2ZwdHJdXG4gIHZhciBmYSA9IGFicyhmaSlcbiAgdmFyIGEsIGJcbiAgaWYoZWEgPCBmYSkge1xuICAgIGIgPSBlaVxuICAgIGVwdHIgKz0gMVxuICAgIGlmKGVwdHIgPCBuZSkge1xuICAgICAgZWkgPSBlW2VwdHJdXG4gICAgICBlYSA9IGFicyhlaSlcbiAgICB9XG4gIH0gZWxzZSB7XG4gICAgYiA9IGZpXG4gICAgZnB0ciArPSAxXG4gICAgaWYoZnB0ciA8IG5mKSB7XG4gICAgICBmaSA9IGZbZnB0cl1cbiAgICAgIGZhID0gYWJzKGZpKVxuICAgIH1cbiAgfVxuICBpZigoZXB0ciA8IG5lICYmIGVhIDwgZmEpIHx8IChmcHRyID49IG5mKSkge1xuICAgIGEgPSBlaVxuICAgIGVwdHIgKz0gMVxuICAgIGlmKGVwdHIgPCBuZSkge1xuICAgICAgZWkgPSBlW2VwdHJdXG4gICAgICBlYSA9IGFicyhlaSlcbiAgICB9XG4gIH0gZWxzZSB7XG4gICAgYSA9IGZpXG4gICAgZnB0ciArPSAxXG4gICAgaWYoZnB0ciA8IG5mKSB7XG4gICAgICBmaSA9IGZbZnB0cl1cbiAgICAgIGZhID0gYWJzKGZpKVxuICAgIH1cbiAgfVxuICB2YXIgeCA9IGEgKyBiXG4gIHZhciBidiA9IHggLSBhXG4gIHZhciB5ID0gYiAtIGJ2XG4gIHZhciBxMCA9IHlcbiAgdmFyIHExID0geFxuICB2YXIgX3gsIF9idiwgX2F2LCBfYnIsIF9hclxuICB3aGlsZShlcHRyIDwgbmUgJiYgZnB0ciA8IG5mKSB7XG4gICAgaWYoZWEgPCBmYSkge1xuICAgICAgYSA9IGVpXG4gICAgICBlcHRyICs9IDFcbiAgICAgIGlmKGVwdHIgPCBuZSkge1xuICAgICAgICBlaSA9IGVbZXB0cl1cbiAgICAgICAgZWEgPSBhYnMoZWkpXG4gICAgICB9XG4gICAgfSBlbHNlIHtcbiAgICAgIGEgPSBmaVxuICAgICAgZnB0ciArPSAxXG4gICAgICBpZihmcHRyIDwgbmYpIHtcbiAgICAgICAgZmkgPSBmW2ZwdHJdXG4gICAgICAgIGZhID0gYWJzKGZpKVxuICAgICAgfVxuICAgIH1cbiAgICBiID0gcTBcbiAgICB4ID0gYSArIGJcbiAgICBidiA9IHggLSBhXG4gICAgeSA9IGIgLSBidlxuICAgIGlmKHkpIHtcbiAgICAgIGdbY291bnQrK10gPSB5XG4gICAgfVxuICAgIF94ID0gcTEgKyB4XG4gICAgX2J2ID0gX3ggLSBxMVxuICAgIF9hdiA9IF94IC0gX2J2XG4gICAgX2JyID0geCAtIF9idlxuICAgIF9hciA9IHExIC0gX2F2XG4gICAgcTAgPSBfYXIgKyBfYnJcbiAgICBxMSA9IF94XG4gIH1cbiAgd2hpbGUoZXB0ciA8IG5lKSB7XG4gICAgYSA9IGVpXG4gICAgYiA9IHEwXG4gICAgeCA9IGEgKyBiXG4gICAgYnYgPSB4IC0gYVxuICAgIHkgPSBiIC0gYnZcbiAgICBpZih5KSB7XG4gICAgICBnW2NvdW50KytdID0geVxuICAgIH1cbiAgICBfeCA9IHExICsgeFxuICAgIF9idiA9IF94IC0gcTFcbiAgICBfYXYgPSBfeCAtIF9idlxuICAgIF9iciA9IHggLSBfYnZcbiAgICBfYXIgPSBxMSAtIF9hdlxuICAgIHEwID0gX2FyICsgX2JyXG4gICAgcTEgPSBfeFxuICAgIGVwdHIgKz0gMVxuICAgIGlmKGVwdHIgPCBuZSkge1xuICAgICAgZWkgPSBlW2VwdHJdXG4gICAgfVxuICB9XG4gIHdoaWxlKGZwdHIgPCBuZikge1xuICAgIGEgPSBmaVxuICAgIGIgPSBxMFxuICAgIHggPSBhICsgYlxuICAgIGJ2ID0geCAtIGFcbiAgICB5ID0gYiAtIGJ2XG4gICAgaWYoeSkge1xuICAgICAgZ1tjb3VudCsrXSA9IHlcbiAgICB9IFxuICAgIF94ID0gcTEgKyB4XG4gICAgX2J2ID0gX3ggLSBxMVxuICAgIF9hdiA9IF94IC0gX2J2XG4gICAgX2JyID0geCAtIF9idlxuICAgIF9hciA9IHExIC0gX2F2XG4gICAgcTAgPSBfYXIgKyBfYnJcbiAgICBxMSA9IF94XG4gICAgZnB0ciArPSAxXG4gICAgaWYoZnB0ciA8IG5mKSB7XG4gICAgICBmaSA9IGZbZnB0cl1cbiAgICB9XG4gIH1cbiAgaWYocTApIHtcbiAgICBnW2NvdW50KytdID0gcTBcbiAgfVxuICBpZihxMSkge1xuICAgIGdbY291bnQrK10gPSBxMVxuICB9XG4gIGlmKCFjb3VudCkge1xuICAgIGdbY291bnQrK10gPSAwLjAgIFxuICB9XG4gIGcubGVuZ3RoID0gY291bnRcbiAgcmV0dXJuIGdcbn0iLCJcInVzZSBzdHJpY3RcIlxuXG5tb2R1bGUuZXhwb3J0cyA9IHR3b1Byb2R1Y3RcblxudmFyIFNQTElUVEVSID0gKyhNYXRoLnBvdygyLCAyNykgKyAxLjApXG5cbmZ1bmN0aW9uIHR3b1Byb2R1Y3QoYSwgYiwgcmVzdWx0KSB7XG4gIHZhciB4ID0gYSAqIGJcblxuICB2YXIgYyA9IFNQTElUVEVSICogYVxuICB2YXIgYWJpZyA9IGMgLSBhXG4gIHZhciBhaGkgPSBjIC0gYWJpZ1xuICB2YXIgYWxvID0gYSAtIGFoaVxuXG4gIHZhciBkID0gU1BMSVRURVIgKiBiXG4gIHZhciBiYmlnID0gZCAtIGJcbiAgdmFyIGJoaSA9IGQgLSBiYmlnXG4gIHZhciBibG8gPSBiIC0gYmhpXG5cbiAgdmFyIGVycjEgPSB4IC0gKGFoaSAqIGJoaSlcbiAgdmFyIGVycjIgPSBlcnIxIC0gKGFsbyAqIGJoaSlcbiAgdmFyIGVycjMgPSBlcnIyIC0gKGFoaSAqIGJsbylcblxuICB2YXIgeSA9IGFsbyAqIGJsbyAtIGVycjNcblxuICBpZihyZXN1bHQpIHtcbiAgICByZXN1bHRbMF0gPSB5XG4gICAgcmVzdWx0WzFdID0geFxuICAgIHJldHVybiByZXN1bHRcbiAgfVxuXG4gIHJldHVybiBbIHksIHggXVxufSIsIlwidXNlIHN0cmljdFwiXG5cbnZhciB0d29Qcm9kdWN0ID0gcmVxdWlyZShcInR3by1wcm9kdWN0XCIpXG52YXIgcm9idXN0U3VtID0gcmVxdWlyZShcInJvYnVzdC1zdW1cIilcbnZhciByb2J1c3RTY2FsZSA9IHJlcXVpcmUoXCJyb2J1c3Qtc2NhbGVcIilcbnZhciByb2J1c3RTdWJ0cmFjdCA9IHJlcXVpcmUoXCJyb2J1c3Qtc3VidHJhY3RcIilcblxudmFyIE5VTV9FWFBBTkQgPSA1XG5cbnZhciBFUFNJTE9OICAgICA9IDEuMTEwMjIzMDI0NjI1MTU2NWUtMTZcbnZhciBFUlJCT1VORDMgICA9ICgzLjAgKyAxNi4wICogRVBTSUxPTikgKiBFUFNJTE9OXG52YXIgRVJSQk9VTkQ0ICAgPSAoNy4wICsgNTYuMCAqIEVQU0lMT04pICogRVBTSUxPTlxuXG5mdW5jdGlvbiBjb2ZhY3RvcihtLCBjKSB7XG4gIHZhciByZXN1bHQgPSBuZXcgQXJyYXkobS5sZW5ndGgtMSlcbiAgZm9yKHZhciBpPTE7IGk8bS5sZW5ndGg7ICsraSkge1xuICAgIHZhciByID0gcmVzdWx0W2ktMV0gPSBuZXcgQXJyYXkobS5sZW5ndGgtMSlcbiAgICBmb3IodmFyIGo9MCxrPTA7IGo8bS5sZW5ndGg7ICsraikge1xuICAgICAgaWYoaiA9PT0gYykge1xuICAgICAgICBjb250aW51ZVxuICAgICAgfVxuICAgICAgcltrKytdID0gbVtpXVtqXVxuICAgIH1cbiAgfVxuICByZXR1cm4gcmVzdWx0XG59XG5cbmZ1bmN0aW9uIG1hdHJpeChuKSB7XG4gIHZhciByZXN1bHQgPSBuZXcgQXJyYXkobilcbiAgZm9yKHZhciBpPTA7IGk8bjsgKytpKSB7XG4gICAgcmVzdWx0W2ldID0gbmV3IEFycmF5KG4pXG4gICAgZm9yKHZhciBqPTA7IGo8bjsgKytqKSB7XG4gICAgICByZXN1bHRbaV1bal0gPSBbXCJtXCIsIGosIFwiW1wiLCAobi1pLTEpLCBcIl1cIl0uam9pbihcIlwiKVxuICAgIH1cbiAgfVxuICByZXR1cm4gcmVzdWx0XG59XG5cbmZ1bmN0aW9uIHNpZ24obikge1xuICBpZihuICYgMSkge1xuICAgIHJldHVybiBcIi1cIlxuICB9XG4gIHJldHVybiBcIlwiXG59XG5cbmZ1bmN0aW9uIGdlbmVyYXRlU3VtKGV4cHIpIHtcbiAgaWYoZXhwci5sZW5ndGggPT09IDEpIHtcbiAgICByZXR1cm4gZXhwclswXVxuICB9IGVsc2UgaWYoZXhwci5sZW5ndGggPT09IDIpIHtcbiAgICByZXR1cm4gW1wic3VtKFwiLCBleHByWzBdLCBcIixcIiwgZXhwclsxXSwgXCIpXCJdLmpvaW4oXCJcIilcbiAgfSBlbHNlIHtcbiAgICB2YXIgbSA9IGV4cHIubGVuZ3RoPj4xXG4gICAgcmV0dXJuIFtcInN1bShcIiwgZ2VuZXJhdGVTdW0oZXhwci5zbGljZSgwLCBtKSksIFwiLFwiLCBnZW5lcmF0ZVN1bShleHByLnNsaWNlKG0pKSwgXCIpXCJdLmpvaW4oXCJcIilcbiAgfVxufVxuXG5mdW5jdGlvbiBkZXRlcm1pbmFudChtKSB7XG4gIGlmKG0ubGVuZ3RoID09PSAyKSB7XG4gICAgcmV0dXJuIFtbXCJzdW0ocHJvZChcIiwgbVswXVswXSwgXCIsXCIsIG1bMV1bMV0sIFwiKSxwcm9kKC1cIiwgbVswXVsxXSwgXCIsXCIsIG1bMV1bMF0sIFwiKSlcIl0uam9pbihcIlwiKV1cbiAgfSBlbHNlIHtcbiAgICB2YXIgZXhwciA9IFtdXG4gICAgZm9yKHZhciBpPTA7IGk8bS5sZW5ndGg7ICsraSkge1xuICAgICAgZXhwci5wdXNoKFtcInNjYWxlKFwiLCBnZW5lcmF0ZVN1bShkZXRlcm1pbmFudChjb2ZhY3RvcihtLCBpKSkpLCBcIixcIiwgc2lnbihpKSwgbVswXVtpXSwgXCIpXCJdLmpvaW4oXCJcIikpXG4gICAgfVxuICAgIHJldHVybiBleHByXG4gIH1cbn1cblxuZnVuY3Rpb24gb3JpZW50YXRpb24obikge1xuICB2YXIgcG9zID0gW11cbiAgdmFyIG5lZyA9IFtdXG4gIHZhciBtID0gbWF0cml4KG4pXG4gIHZhciBhcmdzID0gW11cbiAgZm9yKHZhciBpPTA7IGk8bjsgKytpKSB7XG4gICAgaWYoKGkmMSk9PT0wKSB7XG4gICAgICBwb3MucHVzaC5hcHBseShwb3MsIGRldGVybWluYW50KGNvZmFjdG9yKG0sIGkpKSlcbiAgICB9IGVsc2Uge1xuICAgICAgbmVnLnB1c2guYXBwbHkobmVnLCBkZXRlcm1pbmFudChjb2ZhY3RvcihtLCBpKSkpXG4gICAgfVxuICAgIGFyZ3MucHVzaChcIm1cIiArIGkpXG4gIH1cbiAgdmFyIHBvc0V4cHIgPSBnZW5lcmF0ZVN1bShwb3MpXG4gIHZhciBuZWdFeHByID0gZ2VuZXJhdGVTdW0obmVnKVxuICB2YXIgZnVuY05hbWUgPSBcIm9yaWVudGF0aW9uXCIgKyBuICsgXCJFeGFjdFwiXG4gIHZhciBjb2RlID0gW1wiZnVuY3Rpb24gXCIsIGZ1bmNOYW1lLCBcIihcIiwgYXJncy5qb2luKCksIFwiKXt2YXIgcD1cIiwgcG9zRXhwciwgXCIsbj1cIiwgbmVnRXhwciwgXCIsZD1zdWIocCxuKTtcXFxucmV0dXJuIGRbZC5sZW5ndGgtMV07fTtyZXR1cm4gXCIsIGZ1bmNOYW1lXS5qb2luKFwiXCIpXG4gIHZhciBwcm9jID0gbmV3IEZ1bmN0aW9uKFwic3VtXCIsIFwicHJvZFwiLCBcInNjYWxlXCIsIFwic3ViXCIsIGNvZGUpXG4gIHJldHVybiBwcm9jKHJvYnVzdFN1bSwgdHdvUHJvZHVjdCwgcm9idXN0U2NhbGUsIHJvYnVzdFN1YnRyYWN0KVxufVxuXG52YXIgb3JpZW50YXRpb24zRXhhY3QgPSBvcmllbnRhdGlvbigzKVxudmFyIG9yaWVudGF0aW9uNEV4YWN0ID0gb3JpZW50YXRpb24oNClcblxudmFyIENBQ0hFRCA9IFtcbiAgZnVuY3Rpb24gb3JpZW50YXRpb24wKCkgeyByZXR1cm4gMCB9LFxuICBmdW5jdGlvbiBvcmllbnRhdGlvbjEoKSB7IHJldHVybiAwIH0sXG4gIGZ1bmN0aW9uIG9yaWVudGF0aW9uMihhLCBiKSB7IFxuICAgIHJldHVybiBiWzBdIC0gYVswXVxuICB9LFxuICBmdW5jdGlvbiBvcmllbnRhdGlvbjMoYSwgYiwgYykge1xuICAgIHZhciBsID0gKGFbMV0gLSBjWzFdKSAqIChiWzBdIC0gY1swXSlcbiAgICB2YXIgciA9IChhWzBdIC0gY1swXSkgKiAoYlsxXSAtIGNbMV0pXG4gICAgdmFyIGRldCA9IGwgLSByXG4gICAgdmFyIHNcbiAgICBpZihsID4gMCkge1xuICAgICAgaWYociA8PSAwKSB7XG4gICAgICAgIHJldHVybiBkZXRcbiAgICAgIH0gZWxzZSB7XG4gICAgICAgIHMgPSBsICsgclxuICAgICAgfVxuICAgIH0gZWxzZSBpZihsIDwgMCkge1xuICAgICAgaWYociA+PSAwKSB7XG4gICAgICAgIHJldHVybiBkZXRcbiAgICAgIH0gZWxzZSB7XG4gICAgICAgIHMgPSAtKGwgKyByKVxuICAgICAgfVxuICAgIH0gZWxzZSB7XG4gICAgICByZXR1cm4gZGV0XG4gICAgfVxuICAgIHZhciB0b2wgPSBFUlJCT1VORDMgKiBzXG4gICAgaWYoZGV0ID49IHRvbCB8fCBkZXQgPD0gLXRvbCkge1xuICAgICAgcmV0dXJuIGRldFxuICAgIH1cbiAgICByZXR1cm4gb3JpZW50YXRpb24zRXhhY3QoYSwgYiwgYylcbiAgfSxcbiAgZnVuY3Rpb24gb3JpZW50YXRpb240KGEsYixjLGQpIHtcbiAgICB2YXIgYWR4ID0gYVswXSAtIGRbMF1cbiAgICB2YXIgYmR4ID0gYlswXSAtIGRbMF1cbiAgICB2YXIgY2R4ID0gY1swXSAtIGRbMF1cbiAgICB2YXIgYWR5ID0gYVsxXSAtIGRbMV1cbiAgICB2YXIgYmR5ID0gYlsxXSAtIGRbMV1cbiAgICB2YXIgY2R5ID0gY1sxXSAtIGRbMV1cbiAgICB2YXIgYWR6ID0gYVsyXSAtIGRbMl1cbiAgICB2YXIgYmR6ID0gYlsyXSAtIGRbMl1cbiAgICB2YXIgY2R6ID0gY1syXSAtIGRbMl1cbiAgICB2YXIgYmR4Y2R5ID0gYmR4ICogY2R5XG4gICAgdmFyIGNkeGJkeSA9IGNkeCAqIGJkeVxuICAgIHZhciBjZHhhZHkgPSBjZHggKiBhZHlcbiAgICB2YXIgYWR4Y2R5ID0gYWR4ICogY2R5XG4gICAgdmFyIGFkeGJkeSA9IGFkeCAqIGJkeVxuICAgIHZhciBiZHhhZHkgPSBiZHggKiBhZHlcbiAgICB2YXIgZGV0ID0gYWR6ICogKGJkeGNkeSAtIGNkeGJkeSkgXG4gICAgICAgICAgICArIGJkeiAqIChjZHhhZHkgLSBhZHhjZHkpXG4gICAgICAgICAgICArIGNkeiAqIChhZHhiZHkgLSBiZHhhZHkpXG4gICAgdmFyIHBlcm1hbmVudCA9IChNYXRoLmFicyhiZHhjZHkpICsgTWF0aC5hYnMoY2R4YmR5KSkgKiBNYXRoLmFicyhhZHopXG4gICAgICAgICAgICAgICAgICArIChNYXRoLmFicyhjZHhhZHkpICsgTWF0aC5hYnMoYWR4Y2R5KSkgKiBNYXRoLmFicyhiZHopXG4gICAgICAgICAgICAgICAgICArIChNYXRoLmFicyhhZHhiZHkpICsgTWF0aC5hYnMoYmR4YWR5KSkgKiBNYXRoLmFicyhjZHopXG4gICAgdmFyIHRvbCA9IEVSUkJPVU5ENCAqIHBlcm1hbmVudFxuICAgIGlmICgoZGV0ID4gdG9sKSB8fCAoLWRldCA+IHRvbCkpIHtcbiAgICAgIHJldHVybiBkZXRcbiAgICB9XG4gICAgcmV0dXJuIG9yaWVudGF0aW9uNEV4YWN0KGEsYixjLGQpXG4gIH1cbl1cblxuZnVuY3Rpb24gc2xvd09yaWVudChhcmdzKSB7XG4gIHZhciBwcm9jID0gQ0FDSEVEW2FyZ3MubGVuZ3RoXVxuICBpZighcHJvYykge1xuICAgIHByb2MgPSBDQUNIRURbYXJncy5sZW5ndGhdID0gb3JpZW50YXRpb24oYXJncy5sZW5ndGgpXG4gIH1cbiAgcmV0dXJuIHByb2MuYXBwbHkodW5kZWZpbmVkLCBhcmdzKVxufVxuXG5mdW5jdGlvbiBnZW5lcmF0ZU9yaWVudGF0aW9uUHJvYygpIHtcbiAgd2hpbGUoQ0FDSEVELmxlbmd0aCA8PSBOVU1fRVhQQU5EKSB7XG4gICAgQ0FDSEVELnB1c2gob3JpZW50YXRpb24oQ0FDSEVELmxlbmd0aCkpXG4gIH1cbiAgdmFyIGFyZ3MgPSBbXVxuICB2YXIgcHJvY0FyZ3MgPSBbXCJzbG93XCJdXG4gIGZvcih2YXIgaT0wOyBpPD1OVU1fRVhQQU5EOyArK2kpIHtcbiAgICBhcmdzLnB1c2goXCJhXCIgKyBpKVxuICAgIHByb2NBcmdzLnB1c2goXCJvXCIgKyBpKVxuICB9XG4gIHZhciBjb2RlID0gW1xuICAgIFwiZnVuY3Rpb24gZ2V0T3JpZW50YXRpb24oXCIsIGFyZ3Muam9pbigpLCBcIil7c3dpdGNoKGFyZ3VtZW50cy5sZW5ndGgpe2Nhc2UgMDpjYXNlIDE6cmV0dXJuIDA7XCJcbiAgXVxuICBmb3IodmFyIGk9MjsgaTw9TlVNX0VYUEFORDsgKytpKSB7XG4gICAgY29kZS5wdXNoKFwiY2FzZSBcIiwgaSwgXCI6cmV0dXJuIG9cIiwgaSwgXCIoXCIsIGFyZ3Muc2xpY2UoMCwgaSkuam9pbigpLCBcIik7XCIpXG4gIH1cbiAgY29kZS5wdXNoKFwifXZhciBzPW5ldyBBcnJheShhcmd1bWVudHMubGVuZ3RoKTtmb3IodmFyIGk9MDtpPGFyZ3VtZW50cy5sZW5ndGg7KytpKXtzW2ldPWFyZ3VtZW50c1tpXX07cmV0dXJuIHNsb3cocyk7fXJldHVybiBnZXRPcmllbnRhdGlvblwiKVxuICBwcm9jQXJncy5wdXNoKGNvZGUuam9pbihcIlwiKSlcblxuICB2YXIgcHJvYyA9IEZ1bmN0aW9uLmFwcGx5KHVuZGVmaW5lZCwgcHJvY0FyZ3MpXG4gIG1vZHVsZS5leHBvcnRzID0gcHJvYy5hcHBseSh1bmRlZmluZWQsIFtzbG93T3JpZW50XS5jb25jYXQoQ0FDSEVEKSlcbiAgZm9yKHZhciBpPTA7IGk8PU5VTV9FWFBBTkQ7ICsraSkge1xuICAgIG1vZHVsZS5leHBvcnRzW2ldID0gQ0FDSEVEW2ldXG4gIH1cbn1cblxuZ2VuZXJhdGVPcmllbnRhdGlvblByb2MoKSIsIi8qKlxuICogQml0IHR3aWRkbGluZyBoYWNrcyBmb3IgSmF2YVNjcmlwdC5cbiAqXG4gKiBBdXRob3I6IE1pa29sYSBMeXNlbmtvXG4gKlxuICogUG9ydGVkIGZyb20gU3RhbmZvcmQgYml0IHR3aWRkbGluZyBoYWNrIGxpYnJhcnk6XG4gKiAgICBodHRwOi8vZ3JhcGhpY3Muc3RhbmZvcmQuZWR1L35zZWFuZGVyL2JpdGhhY2tzLmh0bWxcbiAqL1xuXG5cInVzZSBzdHJpY3RcIjsgXCJ1c2UgcmVzdHJpY3RcIjtcblxuLy9OdW1iZXIgb2YgYml0cyBpbiBhbiBpbnRlZ2VyXG52YXIgSU5UX0JJVFMgPSAzMjtcblxuLy9Db25zdGFudHNcbmV4cG9ydHMuSU5UX0JJVFMgID0gSU5UX0JJVFM7XG5leHBvcnRzLklOVF9NQVggICA9ICAweDdmZmZmZmZmO1xuZXhwb3J0cy5JTlRfTUlOICAgPSAtMTw8KElOVF9CSVRTLTEpO1xuXG4vL1JldHVybnMgLTEsIDAsICsxIGRlcGVuZGluZyBvbiBzaWduIG9mIHhcbmV4cG9ydHMuc2lnbiA9IGZ1bmN0aW9uKHYpIHtcbiAgcmV0dXJuICh2ID4gMCkgLSAodiA8IDApO1xufVxuXG4vL0NvbXB1dGVzIGFic29sdXRlIHZhbHVlIG9mIGludGVnZXJcbmV4cG9ydHMuYWJzID0gZnVuY3Rpb24odikge1xuICB2YXIgbWFzayA9IHYgPj4gKElOVF9CSVRTLTEpO1xuICByZXR1cm4gKHYgXiBtYXNrKSAtIG1hc2s7XG59XG5cbi8vQ29tcHV0ZXMgbWluaW11bSBvZiBpbnRlZ2VycyB4IGFuZCB5XG5leHBvcnRzLm1pbiA9IGZ1bmN0aW9uKHgsIHkpIHtcbiAgcmV0dXJuIHkgXiAoKHggXiB5KSAmIC0oeCA8IHkpKTtcbn1cblxuLy9Db21wdXRlcyBtYXhpbXVtIG9mIGludGVnZXJzIHggYW5kIHlcbmV4cG9ydHMubWF4ID0gZnVuY3Rpb24oeCwgeSkge1xuICByZXR1cm4geCBeICgoeCBeIHkpICYgLSh4IDwgeSkpO1xufVxuXG4vL0NoZWNrcyBpZiBhIG51bWJlciBpcyBhIHBvd2VyIG9mIHR3b1xuZXhwb3J0cy5pc1BvdzIgPSBmdW5jdGlvbih2KSB7XG4gIHJldHVybiAhKHYgJiAodi0xKSkgJiYgKCEhdik7XG59XG5cbi8vQ29tcHV0ZXMgbG9nIGJhc2UgMiBvZiB2XG5leHBvcnRzLmxvZzIgPSBmdW5jdGlvbih2KSB7XG4gIHZhciByLCBzaGlmdDtcbiAgciA9ICAgICAodiA+IDB4RkZGRikgPDwgNDsgdiA+Pj49IHI7XG4gIHNoaWZ0ID0gKHYgPiAweEZGICApIDw8IDM7IHYgPj4+PSBzaGlmdDsgciB8PSBzaGlmdDtcbiAgc2hpZnQgPSAodiA+IDB4RiAgICkgPDwgMjsgdiA+Pj49IHNoaWZ0OyByIHw9IHNoaWZ0O1xuICBzaGlmdCA9ICh2ID4gMHgzICAgKSA8PCAxOyB2ID4+Pj0gc2hpZnQ7IHIgfD0gc2hpZnQ7XG4gIHJldHVybiByIHwgKHYgPj4gMSk7XG59XG5cbi8vQ29tcHV0ZXMgbG9nIGJhc2UgMTAgb2YgdlxuZXhwb3J0cy5sb2cxMCA9IGZ1bmN0aW9uKHYpIHtcbiAgcmV0dXJuICAodiA+PSAxMDAwMDAwMDAwKSA/IDkgOiAodiA+PSAxMDAwMDAwMDApID8gOCA6ICh2ID49IDEwMDAwMDAwKSA/IDcgOlxuICAgICAgICAgICh2ID49IDEwMDAwMDApID8gNiA6ICh2ID49IDEwMDAwMCkgPyA1IDogKHYgPj0gMTAwMDApID8gNCA6XG4gICAgICAgICAgKHYgPj0gMTAwMCkgPyAzIDogKHYgPj0gMTAwKSA/IDIgOiAodiA+PSAxMCkgPyAxIDogMDtcbn1cblxuLy9Db3VudHMgbnVtYmVyIG9mIGJpdHNcbmV4cG9ydHMucG9wQ291bnQgPSBmdW5jdGlvbih2KSB7XG4gIHYgPSB2IC0gKCh2ID4+PiAxKSAmIDB4NTU1NTU1NTUpO1xuICB2ID0gKHYgJiAweDMzMzMzMzMzKSArICgodiA+Pj4gMikgJiAweDMzMzMzMzMzKTtcbiAgcmV0dXJuICgodiArICh2ID4+PiA0KSAmIDB4RjBGMEYwRikgKiAweDEwMTAxMDEpID4+PiAyNDtcbn1cblxuLy9Db3VudHMgbnVtYmVyIG9mIHRyYWlsaW5nIHplcm9zXG5mdW5jdGlvbiBjb3VudFRyYWlsaW5nWmVyb3Modikge1xuICB2YXIgYyA9IDMyO1xuICB2ICY9IC12O1xuICBpZiAodikgYy0tO1xuICBpZiAodiAmIDB4MDAwMEZGRkYpIGMgLT0gMTY7XG4gIGlmICh2ICYgMHgwMEZGMDBGRikgYyAtPSA4O1xuICBpZiAodiAmIDB4MEYwRjBGMEYpIGMgLT0gNDtcbiAgaWYgKHYgJiAweDMzMzMzMzMzKSBjIC09IDI7XG4gIGlmICh2ICYgMHg1NTU1NTU1NSkgYyAtPSAxO1xuICByZXR1cm4gYztcbn1cbmV4cG9ydHMuY291bnRUcmFpbGluZ1plcm9zID0gY291bnRUcmFpbGluZ1plcm9zO1xuXG4vL1JvdW5kcyB0byBuZXh0IHBvd2VyIG9mIDJcbmV4cG9ydHMubmV4dFBvdzIgPSBmdW5jdGlvbih2KSB7XG4gIHYgKz0gdiA9PT0gMDtcbiAgLS12O1xuICB2IHw9IHYgPj4+IDE7XG4gIHYgfD0gdiA+Pj4gMjtcbiAgdiB8PSB2ID4+PiA0O1xuICB2IHw9IHYgPj4+IDg7XG4gIHYgfD0gdiA+Pj4gMTY7XG4gIHJldHVybiB2ICsgMTtcbn1cblxuLy9Sb3VuZHMgZG93biB0byBwcmV2aW91cyBwb3dlciBvZiAyXG5leHBvcnRzLnByZXZQb3cyID0gZnVuY3Rpb24odikge1xuICB2IHw9IHYgPj4+IDE7XG4gIHYgfD0gdiA+Pj4gMjtcbiAgdiB8PSB2ID4+PiA0O1xuICB2IHw9IHYgPj4+IDg7XG4gIHYgfD0gdiA+Pj4gMTY7XG4gIHJldHVybiB2IC0gKHY+Pj4xKTtcbn1cblxuLy9Db21wdXRlcyBwYXJpdHkgb2Ygd29yZFxuZXhwb3J0cy5wYXJpdHkgPSBmdW5jdGlvbih2KSB7XG4gIHYgXj0gdiA+Pj4gMTY7XG4gIHYgXj0gdiA+Pj4gODtcbiAgdiBePSB2ID4+PiA0O1xuICB2ICY9IDB4ZjtcbiAgcmV0dXJuICgweDY5OTYgPj4+IHYpICYgMTtcbn1cblxudmFyIFJFVkVSU0VfVEFCTEUgPSBuZXcgQXJyYXkoMjU2KTtcblxuKGZ1bmN0aW9uKHRhYikge1xuICBmb3IodmFyIGk9MDsgaTwyNTY7ICsraSkge1xuICAgIHZhciB2ID0gaSwgciA9IGksIHMgPSA3O1xuICAgIGZvciAodiA+Pj49IDE7IHY7IHYgPj4+PSAxKSB7XG4gICAgICByIDw8PSAxO1xuICAgICAgciB8PSB2ICYgMTtcbiAgICAgIC0tcztcbiAgICB9XG4gICAgdGFiW2ldID0gKHIgPDwgcykgJiAweGZmO1xuICB9XG59KShSRVZFUlNFX1RBQkxFKTtcblxuLy9SZXZlcnNlIGJpdHMgaW4gYSAzMiBiaXQgd29yZFxuZXhwb3J0cy5yZXZlcnNlID0gZnVuY3Rpb24odikge1xuICByZXR1cm4gIChSRVZFUlNFX1RBQkxFWyB2ICAgICAgICAgJiAweGZmXSA8PCAyNCkgfFxuICAgICAgICAgIChSRVZFUlNFX1RBQkxFWyh2ID4+PiA4KSAgJiAweGZmXSA8PCAxNikgfFxuICAgICAgICAgIChSRVZFUlNFX1RBQkxFWyh2ID4+PiAxNikgJiAweGZmXSA8PCA4KSAgfFxuICAgICAgICAgICBSRVZFUlNFX1RBQkxFWyh2ID4+PiAyNCkgJiAweGZmXTtcbn1cblxuLy9JbnRlcmxlYXZlIGJpdHMgb2YgMiBjb29yZGluYXRlcyB3aXRoIDE2IGJpdHMuICBVc2VmdWwgZm9yIGZhc3QgcXVhZHRyZWUgY29kZXNcbmV4cG9ydHMuaW50ZXJsZWF2ZTIgPSBmdW5jdGlvbih4LCB5KSB7XG4gIHggJj0gMHhGRkZGO1xuICB4ID0gKHggfCAoeCA8PCA4KSkgJiAweDAwRkYwMEZGO1xuICB4ID0gKHggfCAoeCA8PCA0KSkgJiAweDBGMEYwRjBGO1xuICB4ID0gKHggfCAoeCA8PCAyKSkgJiAweDMzMzMzMzMzO1xuICB4ID0gKHggfCAoeCA8PCAxKSkgJiAweDU1NTU1NTU1O1xuXG4gIHkgJj0gMHhGRkZGO1xuICB5ID0gKHkgfCAoeSA8PCA4KSkgJiAweDAwRkYwMEZGO1xuICB5ID0gKHkgfCAoeSA8PCA0KSkgJiAweDBGMEYwRjBGO1xuICB5ID0gKHkgfCAoeSA8PCAyKSkgJiAweDMzMzMzMzMzO1xuICB5ID0gKHkgfCAoeSA8PCAxKSkgJiAweDU1NTU1NTU1O1xuXG4gIHJldHVybiB4IHwgKHkgPDwgMSk7XG59XG5cbi8vRXh0cmFjdHMgdGhlIG50aCBpbnRlcmxlYXZlZCBjb21wb25lbnRcbmV4cG9ydHMuZGVpbnRlcmxlYXZlMiA9IGZ1bmN0aW9uKHYsIG4pIHtcbiAgdiA9ICh2ID4+PiBuKSAmIDB4NTU1NTU1NTU7XG4gIHYgPSAodiB8ICh2ID4+PiAxKSkgICYgMHgzMzMzMzMzMztcbiAgdiA9ICh2IHwgKHYgPj4+IDIpKSAgJiAweDBGMEYwRjBGO1xuICB2ID0gKHYgfCAodiA+Pj4gNCkpICAmIDB4MDBGRjAwRkY7XG4gIHYgPSAodiB8ICh2ID4+PiAxNikpICYgMHgwMDBGRkZGO1xuICByZXR1cm4gKHYgPDwgMTYpID4+IDE2O1xufVxuXG5cbi8vSW50ZXJsZWF2ZSBiaXRzIG9mIDMgY29vcmRpbmF0ZXMsIGVhY2ggd2l0aCAxMCBiaXRzLiAgVXNlZnVsIGZvciBmYXN0IG9jdHJlZSBjb2Rlc1xuZXhwb3J0cy5pbnRlcmxlYXZlMyA9IGZ1bmN0aW9uKHgsIHksIHopIHtcbiAgeCAmPSAweDNGRjtcbiAgeCAgPSAoeCB8ICh4PDwxNikpICYgNDI3ODE5MDMzNTtcbiAgeCAgPSAoeCB8ICh4PDw4KSkgICYgMjUxNzE5Njk1O1xuICB4ICA9ICh4IHwgKHg8PDQpKSAgJiAzMjcyMzU2MDM1O1xuICB4ICA9ICh4IHwgKHg8PDIpKSAgJiAxMjI3MTMzNTEzO1xuXG4gIHkgJj0gMHgzRkY7XG4gIHkgID0gKHkgfCAoeTw8MTYpKSAmIDQyNzgxOTAzMzU7XG4gIHkgID0gKHkgfCAoeTw8OCkpICAmIDI1MTcxOTY5NTtcbiAgeSAgPSAoeSB8ICh5PDw0KSkgICYgMzI3MjM1NjAzNTtcbiAgeSAgPSAoeSB8ICh5PDwyKSkgICYgMTIyNzEzMzUxMztcbiAgeCB8PSAoeSA8PCAxKTtcbiAgXG4gIHogJj0gMHgzRkY7XG4gIHogID0gKHogfCAoejw8MTYpKSAmIDQyNzgxOTAzMzU7XG4gIHogID0gKHogfCAoejw8OCkpICAmIDI1MTcxOTY5NTtcbiAgeiAgPSAoeiB8ICh6PDw0KSkgICYgMzI3MjM1NjAzNTtcbiAgeiAgPSAoeiB8ICh6PDwyKSkgICYgMTIyNzEzMzUxMztcbiAgXG4gIHJldHVybiB4IHwgKHogPDwgMik7XG59XG5cbi8vRXh0cmFjdHMgbnRoIGludGVybGVhdmVkIGNvbXBvbmVudCBvZiBhIDMtdHVwbGVcbmV4cG9ydHMuZGVpbnRlcmxlYXZlMyA9IGZ1bmN0aW9uKHYsIG4pIHtcbiAgdiA9ICh2ID4+PiBuKSAgICAgICAmIDEyMjcxMzM1MTM7XG4gIHYgPSAodiB8ICh2Pj4+MikpICAgJiAzMjcyMzU2MDM1O1xuICB2ID0gKHYgfCAodj4+PjQpKSAgICYgMjUxNzE5Njk1O1xuICB2ID0gKHYgfCAodj4+PjgpKSAgICYgNDI3ODE5MDMzNTtcbiAgdiA9ICh2IHwgKHY+Pj4xNikpICAmIDB4M0ZGO1xuICByZXR1cm4gKHY8PDIyKT4+MjI7XG59XG5cbi8vQ29tcHV0ZXMgbmV4dCBjb21iaW5hdGlvbiBpbiBjb2xleGljb2dyYXBoaWMgb3JkZXIgKHRoaXMgaXMgbWlzdGFrZW5seSBjYWxsZWQgbmV4dFBlcm11dGF0aW9uIG9uIHRoZSBiaXQgdHdpZGRsaW5nIGhhY2tzIHBhZ2UpXG5leHBvcnRzLm5leHRDb21iaW5hdGlvbiA9IGZ1bmN0aW9uKHYpIHtcbiAgdmFyIHQgPSB2IHwgKHYgLSAxKTtcbiAgcmV0dXJuICh0ICsgMSkgfCAoKCh+dCAmIC1+dCkgLSAxKSA+Pj4gKGNvdW50VHJhaWxpbmdaZXJvcyh2KSArIDEpKTtcbn1cblxuIiwiXCJ1c2Ugc3RyaWN0XCI7IFwidXNlIHJlc3RyaWN0XCI7XG5cbm1vZHVsZS5leHBvcnRzID0gVW5pb25GaW5kO1xuXG5mdW5jdGlvbiBVbmlvbkZpbmQoY291bnQpIHtcbiAgdGhpcy5yb290cyA9IG5ldyBBcnJheShjb3VudCk7XG4gIHRoaXMucmFua3MgPSBuZXcgQXJyYXkoY291bnQpO1xuICBcbiAgZm9yKHZhciBpPTA7IGk8Y291bnQ7ICsraSkge1xuICAgIHRoaXMucm9vdHNbaV0gPSBpO1xuICAgIHRoaXMucmFua3NbaV0gPSAwO1xuICB9XG59XG5cbnZhciBwcm90byA9IFVuaW9uRmluZC5wcm90b3R5cGVcblxuT2JqZWN0LmRlZmluZVByb3BlcnR5KHByb3RvLCBcImxlbmd0aFwiLCB7XG4gIFwiZ2V0XCI6IGZ1bmN0aW9uKCkge1xuICAgIHJldHVybiB0aGlzLnJvb3RzLmxlbmd0aFxuICB9XG59KVxuXG5wcm90by5tYWtlU2V0ID0gZnVuY3Rpb24oKSB7XG4gIHZhciBuID0gdGhpcy5yb290cy5sZW5ndGg7XG4gIHRoaXMucm9vdHMucHVzaChuKTtcbiAgdGhpcy5yYW5rcy5wdXNoKDApO1xuICByZXR1cm4gbjtcbn1cblxucHJvdG8uZmluZCA9IGZ1bmN0aW9uKHgpIHtcbiAgdmFyIHJvb3RzID0gdGhpcy5yb290cztcbiAgd2hpbGUocm9vdHNbeF0gIT09IHgpIHtcbiAgICB2YXIgeSA9IHJvb3RzW3hdO1xuICAgIHJvb3RzW3hdID0gcm9vdHNbeV07XG4gICAgeCA9IHk7XG4gIH1cbiAgcmV0dXJuIHg7XG59XG5cbnByb3RvLmxpbmsgPSBmdW5jdGlvbih4LCB5KSB7XG4gIHZhciB4ciA9IHRoaXMuZmluZCh4KVxuICAgICwgeXIgPSB0aGlzLmZpbmQoeSk7XG4gIGlmKHhyID09PSB5cikge1xuICAgIHJldHVybjtcbiAgfVxuICB2YXIgcmFua3MgPSB0aGlzLnJhbmtzXG4gICAgLCByb290cyA9IHRoaXMucm9vdHNcbiAgICAsIHhkICAgID0gcmFua3NbeHJdXG4gICAgLCB5ZCAgICA9IHJhbmtzW3lyXTtcbiAgaWYoeGQgPCB5ZCkge1xuICAgIHJvb3RzW3hyXSA9IHlyO1xuICB9IGVsc2UgaWYoeWQgPCB4ZCkge1xuICAgIHJvb3RzW3lyXSA9IHhyO1xuICB9IGVsc2Uge1xuICAgIHJvb3RzW3lyXSA9IHhyO1xuICAgICsrcmFua3NbeHJdO1xuICB9XG59IiwiXCJ1c2Ugc3RyaWN0XCI7IFwidXNlIHJlc3RyaWN0XCI7XG5cbnZhciBiaXRzICAgICAgPSByZXF1aXJlKFwiYml0LXR3aWRkbGVcIilcbiAgLCBVbmlvbkZpbmQgPSByZXF1aXJlKFwidW5pb24tZmluZFwiKVxuXG4vL1JldHVybnMgdGhlIGRpbWVuc2lvbiBvZiBhIGNlbGwgY29tcGxleFxuZnVuY3Rpb24gZGltZW5zaW9uKGNlbGxzKSB7XG4gIHZhciBkID0gMFxuICAgICwgbWF4ID0gTWF0aC5tYXhcbiAgZm9yKHZhciBpPTAsIGlsPWNlbGxzLmxlbmd0aDsgaTxpbDsgKytpKSB7XG4gICAgZCA9IG1heChkLCBjZWxsc1tpXS5sZW5ndGgpXG4gIH1cbiAgcmV0dXJuIGQtMVxufVxuZXhwb3J0cy5kaW1lbnNpb24gPSBkaW1lbnNpb25cblxuLy9Db3VudHMgdGhlIG51bWJlciBvZiB2ZXJ0aWNlcyBpbiBmYWNlc1xuZnVuY3Rpb24gY291bnRWZXJ0aWNlcyhjZWxscykge1xuICB2YXIgdmMgPSAtMVxuICAgICwgbWF4ID0gTWF0aC5tYXhcbiAgZm9yKHZhciBpPTAsIGlsPWNlbGxzLmxlbmd0aDsgaTxpbDsgKytpKSB7XG4gICAgdmFyIGMgPSBjZWxsc1tpXVxuICAgIGZvcih2YXIgaj0wLCBqbD1jLmxlbmd0aDsgajxqbDsgKytqKSB7XG4gICAgICB2YyA9IG1heCh2YywgY1tqXSlcbiAgICB9XG4gIH1cbiAgcmV0dXJuIHZjKzFcbn1cbmV4cG9ydHMuY291bnRWZXJ0aWNlcyA9IGNvdW50VmVydGljZXNcblxuLy9SZXR1cm5zIGEgZGVlcCBjb3B5IG9mIGNlbGxzXG5mdW5jdGlvbiBjbG9uZUNlbGxzKGNlbGxzKSB7XG4gIHZhciBuY2VsbHMgPSBuZXcgQXJyYXkoY2VsbHMubGVuZ3RoKVxuICBmb3IodmFyIGk9MCwgaWw9Y2VsbHMubGVuZ3RoOyBpPGlsOyArK2kpIHtcbiAgICBuY2VsbHNbaV0gPSBjZWxsc1tpXS5zbGljZSgwKVxuICB9XG4gIHJldHVybiBuY2VsbHNcbn1cbmV4cG9ydHMuY2xvbmVDZWxscyA9IGNsb25lQ2VsbHNcblxuLy9SYW5rcyBhIHBhaXIgb2YgY2VsbHMgdXAgdG8gcGVybXV0YXRpb25cbmZ1bmN0aW9uIGNvbXBhcmVDZWxscyhhLCBiKSB7XG4gIHZhciBuID0gYS5sZW5ndGhcbiAgICAsIHQgPSBhLmxlbmd0aCAtIGIubGVuZ3RoXG4gICAgLCBtaW4gPSBNYXRoLm1pblxuICBpZih0KSB7XG4gICAgcmV0dXJuIHRcbiAgfVxuICBzd2l0Y2gobikge1xuICAgIGNhc2UgMDpcbiAgICAgIHJldHVybiAwO1xuICAgIGNhc2UgMTpcbiAgICAgIHJldHVybiBhWzBdIC0gYlswXTtcbiAgICBjYXNlIDI6XG4gICAgICB2YXIgZCA9IGFbMF0rYVsxXS1iWzBdLWJbMV1cbiAgICAgIGlmKGQpIHtcbiAgICAgICAgcmV0dXJuIGRcbiAgICAgIH1cbiAgICAgIHJldHVybiBtaW4oYVswXSxhWzFdKSAtIG1pbihiWzBdLGJbMV0pXG4gICAgY2FzZSAzOlxuICAgICAgdmFyIGwxID0gYVswXSthWzFdXG4gICAgICAgICwgbTEgPSBiWzBdK2JbMV1cbiAgICAgIGQgPSBsMSthWzJdIC0gKG0xK2JbMl0pXG4gICAgICBpZihkKSB7XG4gICAgICAgIHJldHVybiBkXG4gICAgICB9XG4gICAgICB2YXIgbDAgPSBtaW4oYVswXSwgYVsxXSlcbiAgICAgICAgLCBtMCA9IG1pbihiWzBdLCBiWzFdKVxuICAgICAgICAsIGQgID0gbWluKGwwLCBhWzJdKSAtIG1pbihtMCwgYlsyXSlcbiAgICAgIGlmKGQpIHtcbiAgICAgICAgcmV0dXJuIGRcbiAgICAgIH1cbiAgICAgIHJldHVybiBtaW4obDArYVsyXSwgbDEpIC0gbWluKG0wK2JbMl0sIG0xKVxuICAgIFxuICAgIC8vVE9ETzogTWF5YmUgb3B0aW1pemUgbj00IGFzIHdlbGw/XG4gICAgXG4gICAgZGVmYXVsdDpcbiAgICAgIHZhciBhcyA9IGEuc2xpY2UoMClcbiAgICAgIGFzLnNvcnQoKVxuICAgICAgdmFyIGJzID0gYi5zbGljZSgwKVxuICAgICAgYnMuc29ydCgpXG4gICAgICBmb3IodmFyIGk9MDsgaTxuOyArK2kpIHtcbiAgICAgICAgdCA9IGFzW2ldIC0gYnNbaV1cbiAgICAgICAgaWYodCkge1xuICAgICAgICAgIHJldHVybiB0XG4gICAgICAgIH1cbiAgICAgIH1cbiAgICAgIHJldHVybiAwXG4gIH1cbn1cbmV4cG9ydHMuY29tcGFyZUNlbGxzID0gY29tcGFyZUNlbGxzXG5cbmZ1bmN0aW9uIGNvbXBhcmVaaXBwZWQoYSwgYikge1xuICByZXR1cm4gY29tcGFyZUNlbGxzKGFbMF0sIGJbMF0pXG59XG5cbi8vUHV0cyBhIGNlbGwgY29tcGxleCBpbnRvIG5vcm1hbCBvcmRlciBmb3IgdGhlIHB1cnBvc2VzIG9mIGZpbmRDZWxsIHF1ZXJpZXNcbmZ1bmN0aW9uIG5vcm1hbGl6ZShjZWxscywgYXR0cikge1xuICBpZihhdHRyKSB7XG4gICAgdmFyIGxlbiA9IGNlbGxzLmxlbmd0aFxuICAgIHZhciB6aXBwZWQgPSBuZXcgQXJyYXkobGVuKVxuICAgIGZvcih2YXIgaT0wOyBpPGxlbjsgKytpKSB7XG4gICAgICB6aXBwZWRbaV0gPSBbY2VsbHNbaV0sIGF0dHJbaV1dXG4gICAgfVxuICAgIHppcHBlZC5zb3J0KGNvbXBhcmVaaXBwZWQpXG4gICAgZm9yKHZhciBpPTA7IGk8bGVuOyArK2kpIHtcbiAgICAgIGNlbGxzW2ldID0gemlwcGVkW2ldWzBdXG4gICAgICBhdHRyW2ldID0gemlwcGVkW2ldWzFdXG4gICAgfVxuICAgIHJldHVybiBjZWxsc1xuICB9IGVsc2Uge1xuICAgIGNlbGxzLnNvcnQoY29tcGFyZUNlbGxzKVxuICAgIHJldHVybiBjZWxsc1xuICB9XG59XG5leHBvcnRzLm5vcm1hbGl6ZSA9IG5vcm1hbGl6ZVxuXG4vL1JlbW92ZXMgYWxsIGR1cGxpY2F0ZSBjZWxscyBpbiB0aGUgY29tcGxleFxuZnVuY3Rpb24gdW5pcXVlKGNlbGxzKSB7XG4gIGlmKGNlbGxzLmxlbmd0aCA9PT0gMCkge1xuICAgIHJldHVybiBbXVxuICB9XG4gIHZhciBwdHIgPSAxXG4gICAgLCBsZW4gPSBjZWxscy5sZW5ndGhcbiAgZm9yKHZhciBpPTE7IGk8bGVuOyArK2kpIHtcbiAgICB2YXIgYSA9IGNlbGxzW2ldXG4gICAgaWYoY29tcGFyZUNlbGxzKGEsIGNlbGxzW2ktMV0pKSB7XG4gICAgICBpZihpID09PSBwdHIpIHtcbiAgICAgICAgcHRyKytcbiAgICAgICAgY29udGludWVcbiAgICAgIH1cbiAgICAgIGNlbGxzW3B0cisrXSA9IGFcbiAgICB9XG4gIH1cbiAgY2VsbHMubGVuZ3RoID0gcHRyXG4gIHJldHVybiBjZWxsc1xufVxuZXhwb3J0cy51bmlxdWUgPSB1bmlxdWU7XG5cbi8vRmluZHMgYSBjZWxsIGluIGEgbm9ybWFsaXplZCBjZWxsIGNvbXBsZXhcbmZ1bmN0aW9uIGZpbmRDZWxsKGNlbGxzLCBjKSB7XG4gIHZhciBsbyA9IDBcbiAgICAsIGhpID0gY2VsbHMubGVuZ3RoLTFcbiAgICAsIHIgID0gLTFcbiAgd2hpbGUgKGxvIDw9IGhpKSB7XG4gICAgdmFyIG1pZCA9IChsbyArIGhpKSA+PiAxXG4gICAgICAsIHMgICA9IGNvbXBhcmVDZWxscyhjZWxsc1ttaWRdLCBjKVxuICAgIGlmKHMgPD0gMCkge1xuICAgICAgaWYocyA9PT0gMCkge1xuICAgICAgICByID0gbWlkXG4gICAgICB9XG4gICAgICBsbyA9IG1pZCArIDFcbiAgICB9IGVsc2UgaWYocyA+IDApIHtcbiAgICAgIGhpID0gbWlkIC0gMVxuICAgIH1cbiAgfVxuICByZXR1cm4gclxufVxuZXhwb3J0cy5maW5kQ2VsbCA9IGZpbmRDZWxsO1xuXG4vL0J1aWxkcyBhbiBpbmRleCBmb3IgYW4gbi1jZWxsLiAgVGhpcyBpcyBtb3JlIGdlbmVyYWwgdGhhbiBkdWFsLCBidXQgbGVzcyBlZmZpY2llbnRcbmZ1bmN0aW9uIGluY2lkZW5jZShmcm9tX2NlbGxzLCB0b19jZWxscykge1xuICB2YXIgaW5kZXggPSBuZXcgQXJyYXkoZnJvbV9jZWxscy5sZW5ndGgpXG4gIGZvcih2YXIgaT0wLCBpbD1pbmRleC5sZW5ndGg7IGk8aWw7ICsraSkge1xuICAgIGluZGV4W2ldID0gW11cbiAgfVxuICB2YXIgYiA9IFtdXG4gIGZvcih2YXIgaT0wLCBuPXRvX2NlbGxzLmxlbmd0aDsgaTxuOyArK2kpIHtcbiAgICB2YXIgYyA9IHRvX2NlbGxzW2ldXG4gICAgdmFyIGNsID0gYy5sZW5ndGhcbiAgICBmb3IodmFyIGs9MSwga249KDE8PGNsKTsgazxrbjsgKytrKSB7XG4gICAgICBiLmxlbmd0aCA9IGJpdHMucG9wQ291bnQoaylcbiAgICAgIHZhciBsID0gMFxuICAgICAgZm9yKHZhciBqPTA7IGo8Y2w7ICsraikge1xuICAgICAgICBpZihrICYgKDE8PGopKSB7XG4gICAgICAgICAgYltsKytdID0gY1tqXVxuICAgICAgICB9XG4gICAgICB9XG4gICAgICB2YXIgaWR4PWZpbmRDZWxsKGZyb21fY2VsbHMsIGIpXG4gICAgICBpZihpZHggPCAwKSB7XG4gICAgICAgIGNvbnRpbnVlXG4gICAgICB9XG4gICAgICB3aGlsZSh0cnVlKSB7XG4gICAgICAgIGluZGV4W2lkeCsrXS5wdXNoKGkpXG4gICAgICAgIGlmKGlkeCA+PSBmcm9tX2NlbGxzLmxlbmd0aCB8fCBjb21wYXJlQ2VsbHMoZnJvbV9jZWxsc1tpZHhdLCBiKSAhPT0gMCkge1xuICAgICAgICAgIGJyZWFrXG4gICAgICAgIH1cbiAgICAgIH1cbiAgICB9XG4gIH1cbiAgcmV0dXJuIGluZGV4XG59XG5leHBvcnRzLmluY2lkZW5jZSA9IGluY2lkZW5jZVxuXG4vL0NvbXB1dGVzIHRoZSBkdWFsIG9mIHRoZSBtZXNoLiAgVGhpcyBpcyBiYXNpY2FsbHkgYW4gb3B0aW1pemVkIHZlcnNpb24gb2YgYnVpbGRJbmRleCBmb3IgdGhlIHNpdHVhdGlvbiB3aGVyZSBmcm9tX2NlbGxzIGlzIGp1c3QgdGhlIGxpc3Qgb2YgdmVydGljZXNcbmZ1bmN0aW9uIGR1YWwoY2VsbHMsIHZlcnRleF9jb3VudCkge1xuICBpZighdmVydGV4X2NvdW50KSB7XG4gICAgcmV0dXJuIGluY2lkZW5jZSh1bmlxdWUoc2tlbGV0b24oY2VsbHMsIDApKSwgY2VsbHMsIDApXG4gIH1cbiAgdmFyIHJlcyA9IG5ldyBBcnJheSh2ZXJ0ZXhfY291bnQpXG4gIGZvcih2YXIgaT0wOyBpPHZlcnRleF9jb3VudDsgKytpKSB7XG4gICAgcmVzW2ldID0gW11cbiAgfVxuICBmb3IodmFyIGk9MCwgbGVuPWNlbGxzLmxlbmd0aDsgaTxsZW47ICsraSkge1xuICAgIHZhciBjID0gY2VsbHNbaV1cbiAgICBmb3IodmFyIGo9MCwgY2w9Yy5sZW5ndGg7IGo8Y2w7ICsraikge1xuICAgICAgcmVzW2Nbal1dLnB1c2goaSlcbiAgICB9XG4gIH1cbiAgcmV0dXJuIHJlc1xufVxuZXhwb3J0cy5kdWFsID0gZHVhbFxuXG4vL0VudW1lcmF0ZXMgYWxsIGNlbGxzIGluIHRoZSBjb21wbGV4XG5mdW5jdGlvbiBleHBsb2RlKGNlbGxzKSB7XG4gIHZhciByZXN1bHQgPSBbXVxuICBmb3IodmFyIGk9MCwgaWw9Y2VsbHMubGVuZ3RoOyBpPGlsOyArK2kpIHtcbiAgICB2YXIgYyA9IGNlbGxzW2ldXG4gICAgICAsIGNsID0gYy5sZW5ndGh8MFxuICAgIGZvcih2YXIgaj0xLCBqbD0oMTw8Y2wpOyBqPGpsOyArK2opIHtcbiAgICAgIHZhciBiID0gW11cbiAgICAgIGZvcih2YXIgaz0wOyBrPGNsOyArK2spIHtcbiAgICAgICAgaWYoKGogPj4+IGspICYgMSkge1xuICAgICAgICAgIGIucHVzaChjW2tdKVxuICAgICAgICB9XG4gICAgICB9XG4gICAgICByZXN1bHQucHVzaChiKVxuICAgIH1cbiAgfVxuICByZXR1cm4gbm9ybWFsaXplKHJlc3VsdClcbn1cbmV4cG9ydHMuZXhwbG9kZSA9IGV4cGxvZGVcblxuLy9FbnVtZXJhdGVzIGFsbCBvZiB0aGUgbi1jZWxscyBvZiBhIGNlbGwgY29tcGxleFxuZnVuY3Rpb24gc2tlbGV0b24oY2VsbHMsIG4pIHtcbiAgaWYobiA8IDApIHtcbiAgICByZXR1cm4gW11cbiAgfVxuICB2YXIgcmVzdWx0ID0gW11cbiAgICAsIGswICAgICA9ICgxPDwobisxKSktMVxuICBmb3IodmFyIGk9MDsgaTxjZWxscy5sZW5ndGg7ICsraSkge1xuICAgIHZhciBjID0gY2VsbHNbaV1cbiAgICBmb3IodmFyIGs9azA7IGs8KDE8PGMubGVuZ3RoKTsgaz1iaXRzLm5leHRDb21iaW5hdGlvbihrKSkge1xuICAgICAgdmFyIGIgPSBuZXcgQXJyYXkobisxKVxuICAgICAgICAsIGwgPSAwXG4gICAgICBmb3IodmFyIGo9MDsgajxjLmxlbmd0aDsgKytqKSB7XG4gICAgICAgIGlmKGsgJiAoMTw8aikpIHtcbiAgICAgICAgICBiW2wrK10gPSBjW2pdXG4gICAgICAgIH1cbiAgICAgIH1cbiAgICAgIHJlc3VsdC5wdXNoKGIpXG4gICAgfVxuICB9XG4gIHJldHVybiBub3JtYWxpemUocmVzdWx0KVxufVxuZXhwb3J0cy5za2VsZXRvbiA9IHNrZWxldG9uO1xuXG4vL0NvbXB1dGVzIHRoZSBib3VuZGFyeSBvZiBhbGwgY2VsbHMsIGRvZXMgbm90IHJlbW92ZSBkdXBsaWNhdGVzXG5mdW5jdGlvbiBib3VuZGFyeShjZWxscykge1xuICB2YXIgcmVzID0gW11cbiAgZm9yKHZhciBpPTAsaWw9Y2VsbHMubGVuZ3RoOyBpPGlsOyArK2kpIHtcbiAgICB2YXIgYyA9IGNlbGxzW2ldXG4gICAgZm9yKHZhciBqPTAsY2w9Yy5sZW5ndGg7IGo8Y2w7ICsraikge1xuICAgICAgdmFyIGIgPSBuZXcgQXJyYXkoYy5sZW5ndGgtMSlcbiAgICAgIGZvcih2YXIgaz0wLCBsPTA7IGs8Y2w7ICsraykge1xuICAgICAgICBpZihrICE9PSBqKSB7XG4gICAgICAgICAgYltsKytdID0gY1trXVxuICAgICAgICB9XG4gICAgICB9XG4gICAgICByZXMucHVzaChiKVxuICAgIH1cbiAgfVxuICByZXR1cm4gbm9ybWFsaXplKHJlcylcbn1cbmV4cG9ydHMuYm91bmRhcnkgPSBib3VuZGFyeTtcblxuLy9Db21wdXRlcyBjb25uZWN0ZWQgY29tcG9uZW50cyBmb3IgYSBkZW5zZSBjZWxsIGNvbXBsZXhcbmZ1bmN0aW9uIGNvbm5lY3RlZENvbXBvbmVudHNfZGVuc2UoY2VsbHMsIHZlcnRleF9jb3VudCkge1xuICB2YXIgbGFiZWxzID0gbmV3IFVuaW9uRmluZCh2ZXJ0ZXhfY291bnQpXG4gIGZvcih2YXIgaT0wOyBpPGNlbGxzLmxlbmd0aDsgKytpKSB7XG4gICAgdmFyIGMgPSBjZWxsc1tpXVxuICAgIGZvcih2YXIgaj0wOyBqPGMubGVuZ3RoOyArK2opIHtcbiAgICAgIGZvcih2YXIgaz1qKzE7IGs8Yy5sZW5ndGg7ICsraykge1xuICAgICAgICBsYWJlbHMubGluayhjW2pdLCBjW2tdKVxuICAgICAgfVxuICAgIH1cbiAgfVxuICB2YXIgY29tcG9uZW50cyA9IFtdXG4gICAgLCBjb21wb25lbnRfbGFiZWxzID0gbGFiZWxzLnJhbmtzXG4gIGZvcih2YXIgaT0wOyBpPGNvbXBvbmVudF9sYWJlbHMubGVuZ3RoOyArK2kpIHtcbiAgICBjb21wb25lbnRfbGFiZWxzW2ldID0gLTFcbiAgfVxuICBmb3IodmFyIGk9MDsgaTxjZWxscy5sZW5ndGg7ICsraSkge1xuICAgIHZhciBsID0gbGFiZWxzLmZpbmQoY2VsbHNbaV1bMF0pXG4gICAgaWYoY29tcG9uZW50X2xhYmVsc1tsXSA8IDApIHtcbiAgICAgIGNvbXBvbmVudF9sYWJlbHNbbF0gPSBjb21wb25lbnRzLmxlbmd0aFxuICAgICAgY29tcG9uZW50cy5wdXNoKFtjZWxsc1tpXS5zbGljZSgwKV0pXG4gICAgfSBlbHNlIHtcbiAgICAgIGNvbXBvbmVudHNbY29tcG9uZW50X2xhYmVsc1tsXV0ucHVzaChjZWxsc1tpXS5zbGljZSgwKSlcbiAgICB9XG4gIH1cbiAgcmV0dXJuIGNvbXBvbmVudHNcbn1cblxuLy9Db21wdXRlcyBjb25uZWN0ZWQgY29tcG9uZW50cyBmb3IgYSBzcGFyc2UgZ3JhcGhcbmZ1bmN0aW9uIGNvbm5lY3RlZENvbXBvbmVudHNfc3BhcnNlKGNlbGxzKSB7XG4gIHZhciB2ZXJ0aWNlcyAgPSB1bmlxdWUobm9ybWFsaXplKHNrZWxldG9uKGNlbGxzLCAwKSkpXG4gICAgLCBsYWJlbHMgICAgPSBuZXcgVW5pb25GaW5kKHZlcnRpY2VzLmxlbmd0aClcbiAgZm9yKHZhciBpPTA7IGk8Y2VsbHMubGVuZ3RoOyArK2kpIHtcbiAgICB2YXIgYyA9IGNlbGxzW2ldXG4gICAgZm9yKHZhciBqPTA7IGo8Yy5sZW5ndGg7ICsraikge1xuICAgICAgdmFyIHZqID0gZmluZENlbGwodmVydGljZXMsIFtjW2pdXSlcbiAgICAgIGZvcih2YXIgaz1qKzE7IGs8Yy5sZW5ndGg7ICsraykge1xuICAgICAgICBsYWJlbHMubGluayh2aiwgZmluZENlbGwodmVydGljZXMsIFtjW2tdXSkpXG4gICAgICB9XG4gICAgfVxuICB9XG4gIHZhciBjb21wb25lbnRzICAgICAgICA9IFtdXG4gICAgLCBjb21wb25lbnRfbGFiZWxzICA9IGxhYmVscy5yYW5rc1xuICBmb3IodmFyIGk9MDsgaTxjb21wb25lbnRfbGFiZWxzLmxlbmd0aDsgKytpKSB7XG4gICAgY29tcG9uZW50X2xhYmVsc1tpXSA9IC0xXG4gIH1cbiAgZm9yKHZhciBpPTA7IGk8Y2VsbHMubGVuZ3RoOyArK2kpIHtcbiAgICB2YXIgbCA9IGxhYmVscy5maW5kKGZpbmRDZWxsKHZlcnRpY2VzLCBbY2VsbHNbaV1bMF1dKSk7XG4gICAgaWYoY29tcG9uZW50X2xhYmVsc1tsXSA8IDApIHtcbiAgICAgIGNvbXBvbmVudF9sYWJlbHNbbF0gPSBjb21wb25lbnRzLmxlbmd0aFxuICAgICAgY29tcG9uZW50cy5wdXNoKFtjZWxsc1tpXS5zbGljZSgwKV0pXG4gICAgfSBlbHNlIHtcbiAgICAgIGNvbXBvbmVudHNbY29tcG9uZW50X2xhYmVsc1tsXV0ucHVzaChjZWxsc1tpXS5zbGljZSgwKSlcbiAgICB9XG4gIH1cbiAgcmV0dXJuIGNvbXBvbmVudHNcbn1cblxuLy9Db21wdXRlcyBjb25uZWN0ZWQgY29tcG9uZW50cyBmb3IgYSBjZWxsIGNvbXBsZXhcbmZ1bmN0aW9uIGNvbm5lY3RlZENvbXBvbmVudHMoY2VsbHMsIHZlcnRleF9jb3VudCkge1xuICBpZih2ZXJ0ZXhfY291bnQpIHtcbiAgICByZXR1cm4gY29ubmVjdGVkQ29tcG9uZW50c19kZW5zZShjZWxscywgdmVydGV4X2NvdW50KVxuICB9XG4gIHJldHVybiBjb25uZWN0ZWRDb21wb25lbnRzX3NwYXJzZShjZWxscylcbn1cbmV4cG9ydHMuY29ubmVjdGVkQ29tcG9uZW50cyA9IGNvbm5lY3RlZENvbXBvbmVudHNcbiIsIlwidXNlIHN0cmljdFwiXG5cbmZ1bmN0aW9uIHVuaXF1ZV9wcmVkKGxpc3QsIGNvbXBhcmUpIHtcbiAgdmFyIHB0ciA9IDFcbiAgICAsIGxlbiA9IGxpc3QubGVuZ3RoXG4gICAgLCBhPWxpc3RbMF0sIGI9bGlzdFswXVxuICBmb3IodmFyIGk9MTsgaTxsZW47ICsraSkge1xuICAgIGIgPSBhXG4gICAgYSA9IGxpc3RbaV1cbiAgICBpZihjb21wYXJlKGEsIGIpKSB7XG4gICAgICBpZihpID09PSBwdHIpIHtcbiAgICAgICAgcHRyKytcbiAgICAgICAgY29udGludWVcbiAgICAgIH1cbiAgICAgIGxpc3RbcHRyKytdID0gYVxuICAgIH1cbiAgfVxuICBsaXN0Lmxlbmd0aCA9IHB0clxuICByZXR1cm4gbGlzdFxufVxuXG5mdW5jdGlvbiB1bmlxdWVfZXEobGlzdCkge1xuICB2YXIgcHRyID0gMVxuICAgICwgbGVuID0gbGlzdC5sZW5ndGhcbiAgICAsIGE9bGlzdFswXSwgYiA9IGxpc3RbMF1cbiAgZm9yKHZhciBpPTE7IGk8bGVuOyArK2ksIGI9YSkge1xuICAgIGIgPSBhXG4gICAgYSA9IGxpc3RbaV1cbiAgICBpZihhICE9PSBiKSB7XG4gICAgICBpZihpID09PSBwdHIpIHtcbiAgICAgICAgcHRyKytcbiAgICAgICAgY29udGludWVcbiAgICAgIH1cbiAgICAgIGxpc3RbcHRyKytdID0gYVxuICAgIH1cbiAgfVxuICBsaXN0Lmxlbmd0aCA9IHB0clxuICByZXR1cm4gbGlzdFxufVxuXG5mdW5jdGlvbiB1bmlxdWUobGlzdCwgY29tcGFyZSwgc29ydGVkKSB7XG4gIGlmKGxpc3QubGVuZ3RoID09PSAwKSB7XG4gICAgcmV0dXJuIGxpc3RcbiAgfVxuICBpZihjb21wYXJlKSB7XG4gICAgaWYoIXNvcnRlZCkge1xuICAgICAgbGlzdC5zb3J0KGNvbXBhcmUpXG4gICAgfVxuICAgIHJldHVybiB1bmlxdWVfcHJlZChsaXN0LCBjb21wYXJlKVxuICB9XG4gIGlmKCFzb3J0ZWQpIHtcbiAgICBsaXN0LnNvcnQoKVxuICB9XG4gIHJldHVybiB1bmlxdWVfZXEobGlzdClcbn1cblxubW9kdWxlLmV4cG9ydHMgPSB1bmlxdWVcbiIsIlwidXNlIHN0cmljdFwiXG5cbnZhciBjaCA9IHJlcXVpcmUoXCJpbmNyZW1lbnRhbC1jb252ZXgtaHVsbFwiKVxudmFyIHVuaXEgPSByZXF1aXJlKFwidW5pcVwiKVxuXG5tb2R1bGUuZXhwb3J0cyA9IHRyaWFuZ3VsYXRlXG5cbmZ1bmN0aW9uIExpZnRlZFBvaW50KHAsIGkpIHtcbiAgdGhpcy5wb2ludCA9IHBcbiAgdGhpcy5pbmRleCA9IGlcbn1cblxuZnVuY3Rpb24gY29tcGFyZUxpZnRlZChhLCBiKSB7XG4gIHZhciBhcCA9IGEucG9pbnRcbiAgdmFyIGJwID0gYi5wb2ludFxuICB2YXIgZCA9IGFwLmxlbmd0aFxuICBmb3IodmFyIGk9MDsgaTxkOyArK2kpIHtcbiAgICB2YXIgcyA9IGJwW2ldIC0gYXBbaV1cbiAgICBpZihzKSB7XG4gICAgICByZXR1cm4gc1xuICAgIH1cbiAgfVxuICByZXR1cm4gMFxufVxuXG5mdW5jdGlvbiB0cmlhbmd1bGF0ZTFEKG4sIHBvaW50cywgaW5jbHVkZVBvaW50QXRJbmZpbml0eSkge1xuICBpZihuID09PSAxKSB7XG4gICAgaWYoaW5jbHVkZVBvaW50QXRJbmZpbml0eSkge1xuICAgICAgcmV0dXJuIFsgWy0xLCAwXSBdXG4gICAgfSBlbHNlIHtcbiAgICAgIHJldHVybiBbXVxuICAgIH1cbiAgfVxuICB2YXIgbGlmdGVkID0gcG9pbnRzLm1hcChmdW5jdGlvbihwLCBpKSB7XG4gICAgcmV0dXJuIFsgcFswXSwgaSBdXG4gIH0pXG4gIGxpZnRlZC5zb3J0KGZ1bmN0aW9uKGEsYikge1xuICAgIHJldHVybiBhWzBdIC0gYlswXVxuICB9KVxuICB2YXIgY2VsbHMgPSBuZXcgQXJyYXkobiAtIDEpXG4gIGZvcih2YXIgaT0xOyBpPG47ICsraSkge1xuICAgIHZhciBhID0gbGlmdGVkW2ktMV1cbiAgICB2YXIgYiA9IGxpZnRlZFtpXVxuICAgIGNlbGxzW2ktMV0gPSBbIGFbMV0sIGJbMV0gXVxuICB9XG4gIGlmKGluY2x1ZGVQb2ludEF0SW5maW5pdHkpIHtcbiAgICBjZWxscy5wdXNoKFxuICAgICAgWyAtMSwgY2VsbHNbMF1bMV0sIF0sXG4gICAgICBbIGNlbGxzW24tMV1bMV0sIC0xIF0pXG4gIH1cbiAgcmV0dXJuIGNlbGxzXG59XG5cbmZ1bmN0aW9uIHRyaWFuZ3VsYXRlKHBvaW50cywgaW5jbHVkZVBvaW50QXRJbmZpbml0eSkge1xuICB2YXIgbiA9IHBvaW50cy5sZW5ndGhcbiAgaWYobiA9PT0gMCkge1xuICAgIHJldHVybiBbXVxuICB9XG4gIFxuICB2YXIgZCA9IHBvaW50c1swXS5sZW5ndGhcbiAgaWYoZCA8IDEpIHtcbiAgICByZXR1cm4gW11cbiAgfVxuXG4gIC8vU3BlY2lhbCBjYXNlOiAgRm9yIDFEIHdlIGNhbiBqdXN0IHNvcnQgdGhlIHBvaW50c1xuICBpZihkID09PSAxKSB7XG4gICAgcmV0dXJuIHRyaWFuZ3VsYXRlMUQobiwgcG9pbnRzLCBpbmNsdWRlUG9pbnRBdEluZmluaXR5KVxuICB9XG4gIFxuICAvL0xpZnQgcG9pbnRzLCBzb3J0XG4gIHZhciBsaWZ0ZWQgPSBuZXcgQXJyYXkobilcbiAgdmFyIHVwcGVyID0gMS4wXG4gIGZvcih2YXIgaT0wOyBpPG47ICsraSkge1xuICAgIHZhciBwID0gcG9pbnRzW2ldXG4gICAgdmFyIHggPSBuZXcgQXJyYXkoZCsxKVxuICAgIHZhciBsID0gMC4wXG4gICAgZm9yKHZhciBqPTA7IGo8ZDsgKytqKSB7XG4gICAgICB2YXIgdiA9IHBbal1cbiAgICAgIHhbal0gPSB2XG4gICAgICBsICs9IHYgKiB2XG4gICAgfVxuICAgIHhbZF0gPSBsXG4gICAgbGlmdGVkW2ldID0gbmV3IExpZnRlZFBvaW50KHgsIGkpXG4gICAgdXBwZXIgPSBNYXRoLm1heChsLCB1cHBlcilcbiAgfVxuICB1bmlxKGxpZnRlZCwgY29tcGFyZUxpZnRlZClcbiAgXG4gIC8vRG91YmxlIHBvaW50c1xuICBuID0gbGlmdGVkLmxlbmd0aFxuXG4gIC8vQ3JlYXRlIG5ldyBsaXN0IG9mIHBvaW50c1xuICB2YXIgZHBvaW50cyA9IG5ldyBBcnJheShuICsgZCArIDEpXG4gIHZhciBkaW5kZXggPSBuZXcgQXJyYXkobiArIGQgKyAxKVxuXG4gIC8vQWRkIHN0ZWluZXIgcG9pbnRzIGF0IHRvcFxuICB2YXIgdSA9IChkKzEpICogKGQrMSkgKiB1cHBlclxuICB2YXIgeSA9IG5ldyBBcnJheShkKzEpXG4gIGZvcih2YXIgaT0wOyBpPD1kOyArK2kpIHtcbiAgICB5W2ldID0gMC4wXG4gIH1cbiAgeVtkXSA9IHVcblxuICBkcG9pbnRzWzBdID0geS5zbGljZSgpXG4gIGRpbmRleFswXSA9IC0xXG5cbiAgZm9yKHZhciBpPTA7IGk8PWQ7ICsraSkge1xuICAgIHZhciB4ID0geS5zbGljZSgpXG4gICAgeFtpXSA9IDFcbiAgICBkcG9pbnRzW2krMV0gPSB4XG4gICAgZGluZGV4W2krMV0gPSAtMVxuICB9XG5cbiAgLy9Db3B5IHJlc3Qgb2YgdGhlIHBvaW50cyBvdmVyXG4gIGZvcih2YXIgaT0wOyBpPG47ICsraSkge1xuICAgIHZhciBoID0gbGlmdGVkW2ldXG4gICAgZHBvaW50c1tpICsgZCArIDFdID0gaC5wb2ludFxuICAgIGRpbmRleFtpICsgZCArIDFdID0gIGguaW5kZXhcbiAgfVxuXG4gIC8vQ29uc3RydWN0IGNvbnZleCBodWxsXG4gIHZhciBodWxsID0gY2goZHBvaW50cywgZmFsc2UpXG4gIGlmKGluY2x1ZGVQb2ludEF0SW5maW5pdHkpIHtcbiAgICBodWxsID0gaHVsbC5maWx0ZXIoZnVuY3Rpb24oY2VsbCkge1xuICAgICAgdmFyIGNvdW50ID0gMFxuICAgICAgZm9yKHZhciBqPTA7IGo8PWQ7ICsraikge1xuICAgICAgICB2YXIgdiA9IGRpbmRleFtjZWxsW2pdXVxuICAgICAgICBpZih2IDwgMCkge1xuICAgICAgICAgIGlmKCsrY291bnQgPj0gMikge1xuICAgICAgICAgICAgcmV0dXJuIGZhbHNlXG4gICAgICAgICAgfVxuICAgICAgICB9XG4gICAgICAgIGNlbGxbal0gPSB2XG4gICAgICB9XG4gICAgICByZXR1cm4gdHJ1ZVxuICAgIH0pXG4gIH0gZWxzZSB7XG4gICAgaHVsbCA9IGh1bGwuZmlsdGVyKGZ1bmN0aW9uKGNlbGwpIHtcbiAgICAgIGZvcih2YXIgaT0wOyBpPD1kOyArK2kpIHtcbiAgICAgICAgdmFyIHYgPSBkaW5kZXhbY2VsbFtpXV1cbiAgICAgICAgaWYodiA8IDApIHtcbiAgICAgICAgICByZXR1cm4gZmFsc2VcbiAgICAgICAgfVxuICAgICAgICBjZWxsW2ldID0gdlxuICAgICAgfVxuICAgICAgcmV0dXJuIHRydWVcbiAgICB9KVxuICB9XG5cbiAgaWYoZCAmIDEpIHtcbiAgICBmb3IodmFyIGk9MDsgaTxodWxsLmxlbmd0aDsgKytpKSB7XG4gICAgICB2YXIgaCA9IGh1bGxbaV1cbiAgICAgIHZhciB4ID0gaFswXVxuICAgICAgaFswXSA9IGhbMV1cbiAgICAgIGhbMV0gPSB4XG4gICAgfVxuICB9XG5cbiAgcmV0dXJuIGh1bGxcbn0iLCJcbm1vZHVsZS5leHBvcnRzID0gcGFyc2VcblxuLyoqXG4gKiBleHBlY3RlZCBhcmd1bWVudCBsZW5ndGhzXG4gKiBAdHlwZSB7T2JqZWN0fVxuICovXG5cbnZhciBsZW5ndGggPSB7YTogNywgYzogNiwgaDogMSwgbDogMiwgbTogMiwgcTogNCwgczogNCwgdDogMiwgdjogMSwgejogMH1cblxuLyoqXG4gKiBzZWdtZW50IHBhdHRlcm5cbiAqIEB0eXBlIHtSZWdFeHB9XG4gKi9cblxudmFyIHNlZ21lbnQgPSAvKFthc3R2enFtaGxjXSkoW15hc3R2enFtaGxjXSopL2lnXG5cbi8qKlxuICogcGFyc2UgYW4gc3ZnIHBhdGggZGF0YSBzdHJpbmcuIEdlbmVyYXRlcyBhbiBBcnJheVxuICogb2YgY29tbWFuZHMgd2hlcmUgZWFjaCBjb21tYW5kIGlzIGFuIEFycmF5IG9mIHRoZVxuICogZm9ybSBgW2NvbW1hbmQsIGFyZzEsIGFyZzIsIC4uLl1gXG4gKlxuICogQHBhcmFtIHtTdHJpbmd9IHBhdGhcbiAqIEByZXR1cm4ge0FycmF5fVxuICovXG5cbmZ1bmN0aW9uIHBhcnNlKHBhdGgpIHtcblx0dmFyIGRhdGEgPSBbXVxuXHRwYXRoLnJlcGxhY2Uoc2VnbWVudCwgZnVuY3Rpb24oXywgY29tbWFuZCwgYXJncyl7XG5cdFx0dmFyIHR5cGUgPSBjb21tYW5kLnRvTG93ZXJDYXNlKClcblx0XHRhcmdzID0gcGFyc2VWYWx1ZXMoYXJncylcblxuXHRcdC8vIG92ZXJsb2FkZWQgbW92ZVRvXG5cdFx0aWYgKHR5cGUgPT0gJ20nICYmIGFyZ3MubGVuZ3RoID4gMikge1xuXHRcdFx0ZGF0YS5wdXNoKFtjb21tYW5kXS5jb25jYXQoYXJncy5zcGxpY2UoMCwgMikpKVxuXHRcdFx0dHlwZSA9ICdsJ1xuXHRcdFx0Y29tbWFuZCA9IGNvbW1hbmQgPT0gJ20nID8gJ2wnIDogJ0wnXG5cdFx0fVxuXG5cdFx0d2hpbGUgKHRydWUpIHtcblx0XHRcdGlmIChhcmdzLmxlbmd0aCA9PSBsZW5ndGhbdHlwZV0pIHtcblx0XHRcdFx0YXJncy51bnNoaWZ0KGNvbW1hbmQpXG5cdFx0XHRcdHJldHVybiBkYXRhLnB1c2goYXJncylcblx0XHRcdH1cblx0XHRcdGlmIChhcmdzLmxlbmd0aCA8IGxlbmd0aFt0eXBlXSkgdGhyb3cgbmV3IEVycm9yKCdtYWxmb3JtZWQgcGF0aCBkYXRhJylcblx0XHRcdGRhdGEucHVzaChbY29tbWFuZF0uY29uY2F0KGFyZ3Muc3BsaWNlKDAsIGxlbmd0aFt0eXBlXSkpKVxuXHRcdH1cblx0fSlcblx0cmV0dXJuIGRhdGFcbn1cblxuZnVuY3Rpb24gcGFyc2VWYWx1ZXMoYXJncyl7XG5cdGFyZ3MgPSBhcmdzLm1hdGNoKC8tP1suMC05XSsoPzplWy0rXT9cXGQrKT8vaWcpXG5cdHJldHVybiBhcmdzID8gYXJncy5tYXAoTnVtYmVyKSA6IFtdXG59XG4iLCIndXNlIHN0cmljdCc7XG5cbnZhciBzaWduID0gcmVxdWlyZSgnLi91dGlsaXRpZXMuanMnKS5zaWduO1xudmFyIGNhbGN1bGF0ZURpc3RhbmNlID0gcmVxdWlyZSgnLi91dGlsaXRpZXMuanMnKS5kaXN0YW5jZTtcblxuLy8gdmFyIHBvaW50cyA9IHJlcXVpcmUoJy4vaW5pdGlhbGl6ZVBvaW50cy5qcycpLnBvaW50cztcbi8vIHZhciBjaXR5U2V0ID0gcmVxdWlyZSgnLi9pbml0aWFsaXplUG9pbnRzLmpzJykuY2l0eVNldDtcbi8vIHZhciB0ZXh0UG9pbnRzSWQgPSByZXF1aXJlKCcuL2luaXRpYWxpemVQb2ludHMuanMnKS50ZXh0UG9pbnRzSWQ7XG4vLyB2YXIgcG9zc2libGVTdGFydFBvaW50c0lkID0gcmVxdWlyZSgnLi9pbml0aWFsaXplUG9pbnRzLmpzJykucG9zc2libGVTdGFydFBvaW50c0lkO1xuXG52YXIgbGl2ZU1vdXNlUG9zaXRpb24gPSByZXF1aXJlKCcuL21vdXNlLmpzJyk7XG5cbnZhciBWZWN0b3IgPSByZXF1aXJlKCcuL3ZlY3Rvci5qcycpO1xuXG52YXIgcmFuZG9tID0gTWF0aC5yYW5kb207XG52YXIgZmxvb3IgPSBNYXRoLmZsb29yO1xuLy8gdmFyIFJFUFVMU0lPTiA9IDAuMDU7XG4vLyB2YXIgUkVQVUxTSU9OU1BFRUQgPSAwLjAwMjtcbi8vIHZhciBBTlRWRUxPQ0lUWSA9IDAuMDAxO1xuXG5tb2R1bGUuZXhwb3J0cyA9IGZ1bmN0aW9uKGNvbnRhaW5lciwgaW5pdFBvaW50cywgb3B0aW9ucyl7XG5cbiAgICBjb25zb2xlLmxvZygnT3B0aW9ucyBhbnQgOicsIG9wdGlvbnMpO1xuICAgIC8vIERlZmluZSB0aG9zZSBwYXJhbWV0ZXJzIGFzIGF0dHJpYnV0ZXMgb2YgQW50IG9iamVjdCA/XG4gICAgdmFyIFJFUFVMU0lPTiA9IG9wdGlvbnMucmVwU2l6ZTtcbiAgICB2YXIgUkVQVUxTSU9OU1BFRUQgPSBvcHRpb25zLnJlcFNwZWVkO1xuICAgIHZhciBBTlRWRUxPQ0lUWSA9IG9wdGlvbnMudmVsb2NpdHk7XG4gICAgdmFyIFdFSUdIVCA9IG9wdGlvbnMud2VpZ2h0O1xuXG4gICAgdmFyIG1vdXNlID0gbGl2ZU1vdXNlUG9zaXRpb24oY29udGFpbmVyKTtcblxuICAgIHZhciBwb2ludHMgPSBpbml0UG9pbnRzLnBvaW50cztcbiAgICB2YXIgY2l0eVNldCA9IGluaXRQb2ludHMuY2l0eVNldDtcbiAgICB2YXIgdGV4dFBvaW50c0lkID0gaW5pdFBvaW50cy50ZXh0UG9pbnRzSWQ7XG4gICAgdmFyIHBvc3NpYmxlU3RhcnRQb2ludHNJZCA9IGluaXRQb2ludHMucG9zc2libGVTdGFydFBvaW50c0lkO1xuXG5cbiAgICBmdW5jdGlvbiBBbnQocG9pbnQpIHtcbiAgICAgICAgdGhpcy54ID0gcG9pbnQueDsgICAgICAgICAgICAgICAgXG4gICAgICAgIHRoaXMueSA9IHBvaW50Lnk7XG4gICAgICAgIHRoaXMudmVsb2NpdHkgPSBBTlRWRUxPQ0lUWTtcbiAgICAgICAgdGhpcy53ZWlnaHQgPSBXRUlHSFQ7XG4gICAgICAgIHRoaXMucmVwU2l6ZSA9IFJFUFVMU0lPTjtcbiAgICAgICAgdGhpcy5yZXBTcGVlZCA9IFJFUFVMU0lPTlNQRUVEO1xuICAgICAgICB0aGlzLmVkZ2UgPSB1bmRlZmluZWQ7XG4gICAgICAgIHRoaXMuc3RhdGUgPSBcImZvcmFnZVwiO1xuICAgICAgICB0aGlzLmVkZ2VzID0gW107XG4gICAgICAgIHRoaXMubGFzdENpdHkgPSB1bmRlZmluZWQ7XG4gICAgICAgIHRoaXMub3JpZ2luID0gcG9pbnQ7XG4gICAgICAgIHRoaXMuZGVzdGluYXRpb24gPSB1bmRlZmluZWQ7XG4gICAgICAgIHRoaXMub3JpZW50YXRpb24gPSB1bmRlZmluZWQ7XG4gICAgICAgIHRoaXMuZGlyZWN0aW9uID0gbmV3IFZlY3RvcigwLDApO1xuICAgICAgICB0aGlzLnByb2cgPSAwO1xuICAgIH1cbiAgICAvLyBmb3JhZ2U6IHRoZSBhbnQgd2FuZGVycyBhcm91bmQgd2l0aG91dCBhbnkgcGhlcm9tb24gZGVwb3NpdGlvblxuICAgIC8vIG9uY2UgaXQgZmluZHMgYSBjaXR5LCBpdCBzdGFydHMgcmVtZW1iZXJpbmcgdGhlIG5vZGVzIGl0IGdvZXMgdGhyb3VnaFxuICAgIC8vIHdoZW4gaXQgZmluZHMgYW5vdGhlciBjaXR5LCBpdCBjb21wdXRlcyB0aGUgcGF0aCBsZW5ndGggYW5kIGFkZHMgcGhlcm9tb25zIG9uZSBlYWNoIGVkZ2VzXG4gICAgLy8gcHJvcG9ydGlvbm5hbHkgdG8gdGhlIHNob3J0ZXN0bmVzcyBvZiB0aGUgcGF0aFxuICAgIC8vIGl0IHJlc2V0cyB0aGUgbGlzdCBvZiBub2RlcyBhbmQgY29udGludWVzXG4gICAgLy8gd2hpbGUgZm9yYWdpbmcgdGhlIGFudCBjaG9zZXMgdGhlIHBhdGggd2l0aCBhIHBoZXJvbW9uIHByZWZlcmVuY2VcblxuXG4gICAgLy8gc3RhdGljIG1ldGhvZHNcbiAgICBBbnQuZ2VuZXJhdGVSYW5kU3RhcnRQb2ludCA9IGZ1bmN0aW9uKCkge1xuICAgICAgICB2YXIgcmFuZElkID0gTWF0aC5mbG9vcihwb3NzaWJsZVN0YXJ0UG9pbnRzSWQubGVuZ3RoICogcmFuZG9tKCkpO1xuICAgICAgICB2YXIgcmFuZFN0YXJ0UG9pbnQgPSBwb2ludHNbcG9zc2libGVTdGFydFBvaW50c0lkW3JhbmRJZF1dO1xuICAgICAgICByZXR1cm4gcmFuZFN0YXJ0UG9pbnQ7XG4gICAgfVxuXG5cbiAgICAvLyBtZXRob2RzXG4gICAgQW50LnByb3RvdHlwZSA9IHtcblxuICAgICAgICB0cmFuc2l0OiBmdW5jdGlvbigpe1xuICAgICAgICAgICAgc3dpdGNoICh0aGlzLnN0YXRlKSB7XG4gICAgICAgICAgICBjYXNlIFwiZm9yYWdlXCI6XG4gICAgICAgICAgICAgICAgdmFyIHJlcyA9IHRoaXMubW92ZSgpO1xuICAgICAgICAgICAgICAgIGlmIChyZXMuY2l0eVJlYWNoZWQpIHtcbiAgICAgICAgICAgICAgICAgICAgdGhpcy5zdGF0ZSA9IFwicGhlcm9tb25pbmdcIjtcbiAgICAgICAgICAgICAgICAgICAgdGhpcy5sYXN0Q2l0eSA9IHRoaXMub3JpZ2luLmlkO1xuICAgICAgICAgICAgICAgIH07XG4gICAgICAgICAgICAgICAgYnJlYWs7XG4gICAgICAgICAgICBjYXNlIFwicGhlcm9tb25pbmdcIjpcbiAgICAgICAgICAgICAgICB2YXIgcmVzID0gdGhpcy5tb3ZlKCk7XG4gICAgICAgICAgICAgICAgaWYgKHJlcy5lZGdlQ2hhbmdlZCkge1xuICAgICAgICAgICAgICAgICAgICB0aGlzLmVkZ2VzLnB1c2godGhpcy5lZGdlKTtcbiAgICAgICAgICAgICAgICAgICAgLy8gZm91bmQgYSBjaXR5XG4gICAgICAgICAgICAgICAgICAgIGlmIChyZXMuY2l0eVJlYWNoZWQgJiYgKHRoaXMub3JpZ2luLmlkICE9IHRoaXMubGFzdENpdHkpICl7XG4gICAgICAgICAgICAgICAgICAgICAgICAvLyBjb21wdXRlIHRoZSBsZW5ndGggb2YgdGhlIHBhdGhcbiAgICAgICAgICAgICAgICAgICAgICAgIHZhciBwYXRoTGVuZ3RoID0gdGhpcy5lZGdlcy5tYXAoZnVuY3Rpb24oZSl7cmV0dXJuIGUuZGlzdGFuY2V9KS5yZWR1Y2UoZnVuY3Rpb24oYSxiKXtyZXR1cm4gYSArIGJ9KTtcbiAgICAgICAgICAgICAgICAgICAgICAgIHZhciBkZWx0YVBoZXJvbW9uZSA9IDEvcGF0aExlbmd0aDtcbiAgICAgICAgICAgICAgICAgICAgICAgIHZhciBhbnRXZWlnaHQgPSB0aGlzLndlaWdodDtcbiAgICAgICAgICAgICAgICAgICAgICAgIHRoaXMuZWRnZXMuZm9yRWFjaChmdW5jdGlvbihlKXtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICB2YXIgYSA9IGUucHQxLCBiID0gZS5wdDIsIHdlaWdodCA9IDE7ICBcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAvLyBpbmNyZWFzZWQgZHJvcHBlZCBwaGVyb21vbnMgZm9yIHRleHRFZGdlc1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgIGlmICgoY2l0eVNldC5pbmRleE9mKGEuaWQpICE9IC0xKSAmJiBjaXR5U2V0LmluZGV4T2YoYi5pZCkgIT0gLTEgJiYgKE1hdGguYWJzKGEuaWQgLSBiLmlkKSA9PSAxKSlcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHdlaWdodCAqPSBhbnRXZWlnaHQ7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgIGUucGhlcm9tb24gKz0gKGRlbHRhUGhlcm9tb25lICogd2VpZ2h0KTtcbiAgICAgICAgICAgICAgICAgICAgICAgIH0pO1xuXG4gICAgICAgICAgICAgICAgICAgICAgICB0aGlzLmVkZ2VzID0gW3RoaXMuZWRnZV07XG4gICAgICAgICAgICAgICAgICAgICAgICB0aGlzLmxhc3RDaXR5ID0gdGhpcy5vcmlnaW4uaWQ7XG4gICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgIGJyZWFrO1xuICAgICAgICAgICAgfVxuXG4gICAgICAgIH0sXG5cbiAgICAgICAgc2V0RGlyZWN0aW9uOiBmdW5jdGlvbigpe1xuICAgICAgICAgICAgdmFyIHBvc3NpYmxlRWRnZXMgPSBbXTtcblxuICAgICAgICAgICAgZm9yICh2YXIgaSA9IDA7IGkgPCB0aGlzLm9yaWdpbi5uZXh0cy5sZW5ndGg7IGkrKylcbiAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICBwb3NzaWJsZUVkZ2VzW2ldID0gdGhpcy5vcmlnaW4ubmV4dHNbaV07XG4gICAgICAgICAgICB9IFxuXG4gICAgICAgICAgICBwb3NzaWJsZUVkZ2VzLnNwbGljZShwb3NzaWJsZUVkZ2VzLmluZGV4T2YodGhpcy5lZGdlKSwxKTtcblxuICAgICAgICAgICAgLy8gZmxpcCBhIGNvaW4gYW5kIGVpdGhlciB0YWtlIHRoZSBzbWVsbGllc3QgcGF0aCBvciBhIHJhbmRvbSBvbmVcbiAgICAgICAgICAgIGlmIChyYW5kb20oKSA+IDAuNSl7XG4gICAgICAgICAgICAgICAgdmFyIHNtZWxscyA9IHBvc3NpYmxlRWRnZXMubWFwKGZ1bmN0aW9uKGUpe3JldHVybiBlLnBoZXJvbW9ufSk7XG4gICAgICAgICAgICAgICAgdmFyIGluZGV4ID0gc21lbGxzLmluZGV4T2YoTWF0aC5tYXguYXBwbHkoTWF0aCwgc21lbGxzKSk7XG4gICAgICAgICAgICAgICAgdGhpcy5lZGdlID0gcG9zc2libGVFZGdlc1tpbmRleF07XG4gICAgICAgICAgICB9IFxuICAgICAgICAgICAgZWxzZVxuICAgICAgICAgICAgICAgIHRoaXMuZWRnZSA9IHBvc3NpYmxlRWRnZXNbZmxvb3IocmFuZG9tKCkqcG9zc2libGVFZGdlcy5sZW5ndGgpXTtcblxuICAgICAgICAgICAgLy8gc2V0IHRoZSBkZXN0aW5hdGlvbiBwb2ludCwgYmVpbmcgZWRnZS5wdDEgb3IgZWRnZS5wdDJcbiAgICAgICAgICAgIHRoaXMuZGVzdGluYXRpb24gPSAodGhpcy5vcmlnaW4gPT0gdGhpcy5lZGdlLnB0MSkgPyB0aGlzLmVkZ2UucHQyIDogdGhpcy5lZGdlLnB0MTtcblxuICAgICAgICAgICAgdGhpcy5kaXJlY3Rpb24ueCA9IHRoaXMuZGVzdGluYXRpb24ueCAtIHRoaXMub3JpZ2luLng7IFxuICAgICAgICAgICAgdGhpcy5kaXJlY3Rpb24ueSA9IHRoaXMuZGVzdGluYXRpb24ueSAtIHRoaXMub3JpZ2luLnk7XG5cbiAgICAgICAgICAgIHRoaXMuZGlyZWN0aW9uLm5vcm1hbGl6ZSgpO1xuICAgICAgICB9LFxuXG4gICAgICAgIG1vdmU6IGZ1bmN0aW9uKCl7XG4gICAgICAgICAgICB2YXIgZWRnZUNoYW5nZWQ7XG4gICAgICAgICAgICB2YXIgY2l0eVJlYWNoZWQgPSBmYWxzZTtcblxuICAgICAgICAgICAgdGhpcy5kaXJlY3Rpb24ueCA9IHRoaXMuZGVzdGluYXRpb24ueCAtIHRoaXMueDsgXG4gICAgICAgICAgICB0aGlzLmRpcmVjdGlvbi55ID0gdGhpcy5kZXN0aW5hdGlvbi55IC0gdGhpcy55O1xuICAgICAgICAgICAgdGhpcy5kaXJlY3Rpb24ubm9ybWFsaXplKCk7XG5cbiAgICAgICAgICAgIC8vIG9uIGVkZ2VcbiAgICAgICAgICAgIGlmICgoY2FsY3VsYXRlRGlzdGFuY2UodGhpcywgdGhpcy5kZXN0aW5hdGlvbikgPiB0aGlzLnJlcFNwZWVkKSl7XG5cbiAgICAgICAgICAgICAgICAvLyBhIGRlbHRhIG1vdmVtZW50IHdpbGwgYmUgYXBwbGllZCBpZiBjb2xsaXNpb24gd2l0aCBvYnN0YWNsZSBkZXRlY3RlZFxuICAgICAgICAgICAgICAgIHZhciBkZWx0YSA9IHRoaXMuYXZvaWRPYnN0YWNsZSgpO1xuXG4gICAgICAgICAgICAgICAgdGhpcy54ICs9IHRoaXMudmVsb2NpdHkgKiB0aGlzLmRpcmVjdGlvbi54ICsgZGVsdGEueCAqIHRoaXMucmVwU3BlZWQ7XG4gICAgICAgICAgICAgICAgdGhpcy55ICs9IHRoaXMudmVsb2NpdHkgKiB0aGlzLmRpcmVjdGlvbi55ICsgZGVsdGEueSAqIHRoaXMucmVwU3BlZWQ7XG5cbiAgICAgICAgICAgICAgICB0aGlzLnByb2cgPSB0aGlzLmNhbGN1bGF0ZVByb2dyZXNzaW9uKCk7XG4gICAgICAgICAgICAgICAgXG4gICAgICAgICAgICAgICAgZWRnZUNoYW5nZWQgPSBmYWxzZTtcblxuICAgICAgICAgICAgLy8gb24gdmVydGV4XG4gICAgICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgICAgIHRoaXMuc3RlcCA9IDA7XG4gICAgICAgICAgICAgICAgdGhpcy5wcm9nID0gMDtcbiAgICAgICAgICAgICAgICB0aGlzLm9yaWdpbiA9IHRoaXMuZGVzdGluYXRpb247XG4gICAgICAgICAgICAgICAgdGhpcy54ID0gdGhpcy5vcmlnaW4ueDtcbiAgICAgICAgICAgICAgICB0aGlzLnkgPSB0aGlzLm9yaWdpbi55O1xuXG4gICAgICAgICAgICAgICAgdGhpcy5zZXREaXJlY3Rpb24oKTtcblxuICAgICAgICAgICAgICAgIGNpdHlSZWFjaGVkID0gKGNpdHlTZXQuaW5kZXhPZih0aGlzLm9yaWdpbi5pZCkgIT0gLTEpO1xuICAgICAgICAgICAgICAgIGVkZ2VDaGFuZ2VkID0gdHJ1ZTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgICAgIHJldHVybiB7Y2l0eVJlYWNoZWQ6IGNpdHlSZWFjaGVkLCBlZGdlQ2hhbmdlZDogZWRnZUNoYW5nZWR9O1xuICAgICAgICB9LFxuXG4gICAgICAgIGF2b2lkT2JzdGFjbGU6IGZ1bmN0aW9uKCl7XG4gICAgICAgICAgICB2YXIgZGlzdGFuY2UgPSBjYWxjdWxhdGVEaXN0YW5jZSh0aGlzLCBtb3VzZSk7XG4gICAgICAgIFxuICAgICAgICAgICAgaWYgKGRpc3RhbmNlIDw9IHRoaXMucmVwU2l6ZSkge1xuXG4gICAgICAgICAgICAgICAgcmV0dXJuIHtcbiAgICAgICAgICAgICAgICAgICAgLy8gZGVsdGEgbW92ZW1lbnQgaXMgY29tcG9zZWQgb2YgYSByZXB1bHNpb24gZGVsdGEgYW5kIGEgY2lyY3VsYXIgZGVsdGEgXG4gICAgICAgICAgICAgICAgICAgIHg6ICh0aGlzLnggLSBtb3VzZS54KS9kaXN0YW5jZSArICh0aGlzLnkgLSBtb3VzZS55KS9kaXN0YW5jZSAqIDEsXG4gICAgICAgICAgICAgICAgICAgIHk6ICh0aGlzLnkgLSBtb3VzZS55KS9kaXN0YW5jZSAtICh0aGlzLnggLSBtb3VzZS54KS9kaXN0YW5jZSAqIDFcbiAgICAgICAgICAgICAgICB9O1xuICAgICAgICAgICAgfVxuICAgICAgICAgICAgZWxzZVxuICAgICAgICAgICAgICAgIHJldHVybiB7eDowLCB5OjB9O1xuICAgICAgICB9LFxuXG4gICAgICAgIGNhbGN1bGF0ZVByb2dyZXNzaW9uOiBmdW5jdGlvbigpe1xuICAgICAgICAgICAgdmFyIHYgPSBuZXcgVmVjdG9yKHRoaXMueCAtIHRoaXMub3JpZ2luLngsIHRoaXMueSAtIHRoaXMub3JpZ2luLnkpO1xuICAgICAgICAgICAgdmFyIG5vcm0gPSB2Lm5vcm0oKTtcblxuICAgICAgICAgICAgdmFyIHRoZXRhID0gKHYueCAqIHRoaXMuZWRnZS5kaXJlY3Rpb24ueCArIHYueSAqIHRoaXMuZWRnZS5kaXJlY3Rpb24ueSkgLyBub3JtO1xuICAgICAgICAgICAgdmFyIHByb2cgPSBub3JtICogTWF0aC5hYnModGhldGEpO1xuICAgICAgICAgICAgLy8gcmV0dXJucyBsZW5ndGggb2YgcHJvamVjdGlvbiBvbiBlZGdlXG4gICAgICAgICAgICByZXR1cm4gcHJvZztcbiAgICAgICAgfVxuXG4gICAgfTtcbiAgICByZXR1cm4gQW50O1xufVxuXG4iLCIndXNlIHN0cmljdCdcblxudmFyIGR0ID0gcmVxdWlyZShcImRlbGF1bmF5LXRyaWFuZ3VsYXRlXCIpO1xuXG52YXIgcmFuZ2UgPSByZXF1aXJlKCcuL3V0aWxpdGllcy5qcycpLnJhbmdlO1xuXG4vLyB2YXIgcG9pbnRzID0gcmVxdWlyZSgnLi9pbml0aWFsaXplUG9pbnRzLmpzJykucG9pbnRzO1xudmFyIHRleHRNZXNoID0gcmVxdWlyZSgnLi9pbml0aWFsaXplUG9pbnRzLmpzJykudGV4dE1lc2g7XG52YXIgY2l0eVNldCA9IHJlcXVpcmUoJy4vaW5pdGlhbGl6ZVBvaW50cy5qcycpLmNpdHlTZXQ7XG52YXIgbmJSYW5kb21Qb2ludHMgPSByZXF1aXJlKCcuL2luaXRpYWxpemVQb2ludHMuanMnKS5uYlJhbmRvbVBvaW50cztcbnZhciBmb3JjZWRFZGdlcyA9IHJlcXVpcmUoJy4vaW5pdGlhbGl6ZVBvaW50cy5qcycpLmZvcmNlZEVkZ2VzO1xuXG52YXIgRWRnZSA9IHJlcXVpcmUoJy4vZWRnZS5qcycpO1xuXG5cbm1vZHVsZS5leHBvcnRzID0gZnVuY3Rpb24ocG9pbnRzKXtcbiAgICAvLyB0cmlhbmd1bGF0ZVxuICAgIHZhciBjZWxscyA9IGR0KHBvaW50cy5tYXAoZnVuY3Rpb24ocCl7XG4gICAgICAgIHJldHVybiBbcC54LCBwLnldO1xuICAgIH0pKTtcblxuICAgIHZhciBlZGdlcyA9IFtdO1xuICAgIHZhciBwZXJtdXRhdGlvbnMgPSBbWzAsMV0sIFswLDJdLCBbMSwyXV07XG5cbiAgICAvLyBmb3JjZSB0aGUgZWRnZXMgb2YgdGhlIHRleHQgdG8gYmUgZWRnZXMgb2YgdGhlIGdyYXBoXG4gICAgaWYgKHRleHRNZXNoKSB7XG4gICAgICAgIHJhbmdlKDAsIHBvaW50cy5sZW5ndGggLSBuYlJhbmRvbVBvaW50cykuZm9yRWFjaChmdW5jdGlvbihpZCl7XG4gICAgICAgICAgICB2YXIgZGlyZWN0TGluayA9IGZvcmNlZEVkZ2VzW2lkXTtcbiAgICAgICAgICAgIHZhciB0ZXh0RWRnZSA9IEVkZ2UuY3JlYXRlKHBvaW50c1tpZF0sIHBvaW50c1tkaXJlY3RMaW5rXSk7XG4gICAgICAgICAgICBlZGdlcy5wdXNoKHRleHRFZGdlKTtcbiAgICAgICAgICAgIHBvaW50c1tpZF0ubmV4dHMucHVzaCh0ZXh0RWRnZSk7XG4gICAgICAgIH0pXG4gICAgfVxuXG5cbiAgICBjZWxscy5mb3JFYWNoKGZ1bmN0aW9uKGNlbGwpe1xuICAgICAgIFxuICAgICAgICBmb3IgKHZhciBpID0gMDsgaSA8IDM7ICsraSl7ICAvLyBmb3IgZWFjaCBwb2ludC5pZCBsaXN0ZWQgaW4gY3VycmVudCBjZWxsXG4gICAgICAgICAgICB2YXIgcHQgPSBwb2ludHNbY2VsbFtpXV07XG5cbiAgICAgICAgICAgIGZvciAodmFyIGogPSAxOyBqIDwgMzsgKytqKXsgXG5cbiAgICAgICAgICAgICAgICB2YXIgcHRqID0gcG9pbnRzW2NlbGxbKCBpICsgaiApICUgM11dOyAvLyBwaWNrIG9uZSBvZiB0aGUgb3RoZXIgMiBwb2ludHMgb2YgdGhlIGNlbGxcbiAgICAgICAgICAgICAgICB2YXIgbmV3RWRnZSA9IHVuZGVmaW5lZDtcblxuICAgICAgICAgICAgICAgIC8vIGlmIHB0IGFscmVhZHkgaGFzIG5leHRFZGdlcyAuLi5cbiAgICAgICAgICAgICAgICBpZiAocHQubmV4dHMubGVuZ3RoICE9IDApIHtcbiAgICAgICAgICAgICAgICAgICAgXG4gICAgICAgICAgICAgICAgICAgIC8vIC4uLiBnZXQgdGhlIHBvaW50cyBjb3JyZXNwb25kaW5nIC4uLlxuICAgICAgICAgICAgICAgICAgICB2YXIgdGVtcFBvaW50cyA9IHB0Lm5leHRzLm1hcChmdW5jdGlvbihlKXtcbiAgICAgICAgICAgICAgICAgICAgICAgIHJldHVybiBbZS5wdDEsIGUucHQyXTtcbiAgICAgICAgICAgICAgICAgICAgfSkucmVkdWNlKGZ1bmN0aW9uKGEsIGIpe1xuICAgICAgICAgICAgICAgICAgICAgICAgIHJldHVybiBhLmNvbmNhdChiKTtcbiAgICAgICAgICAgICAgICAgICAgfSk7XG5cbiAgICAgICAgICAgICAgICAgICAgLy8gLi4uIGFuZCBjaGVjayBpZiBwdGogYWxyZWFkeSBpcyBwYXJ0IG9mIHRoZSBleGlzdGluZyBuZXh0RWRnZXMuIElmIG5vdCwgYWRkIHRoZSBlZGdlLlxuICAgICAgICAgICAgICAgICAgICBpZiAodGVtcFBvaW50cy5pbmRleE9mKHB0aikgPT0gLTEpe1xuICAgICAgICAgICAgICAgICAgICAgICAgbmV3RWRnZSA9IEVkZ2UuY3JlYXRlKHB0LCBwdGopO1xuICAgICAgICAgICAgICAgICAgICAgICAgZWRnZXMucHVzaChuZXdFZGdlKTtcbiAgICAgICAgICAgICAgICAgICAgICAgIHB0Lm5leHRzLnB1c2gobmV3RWRnZSk7XG4gICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgZWxzZSB7XG4gICAgICAgICAgICAgICAgICAgIG5ld0VkZ2UgPSBFZGdlLmNyZWF0ZShwdCwgcHRqKTtcbiAgICAgICAgICAgICAgICAgICAgZWRnZXMucHVzaChuZXdFZGdlKTtcbiAgICAgICAgICAgICAgICAgICAgcHQubmV4dHMucHVzaChuZXdFZGdlKTtcbiAgICAgICAgICAgICAgICB9XG5cbiAgICAgICAgICAgICAgICAvLyBhZGQgYWxzbyB0aGUgZWRnZSB0byB0aGUgZWRnZSdzIG90aGVyIHBvaW50J3MgbmV4dEVkZ2VzXG4gICAgICAgICAgICAgICAgaWYgKG5ld0VkZ2UgIT0gdW5kZWZpbmVkKXtcbiAgICAgICAgICAgICAgICAgICAgcHRqLm5leHRzLnB1c2gobmV3RWRnZSk7XG4gICAgICAgICAgICAgICAgfSAgICAgICAgIFxuICAgICAgICAgICAgfVxuXG4gICAgICAgICAgICAvLyBhZGQgdGhlIHRleHRFZGdlcyB0byBuZXh0RWRnZXMgbWFwXG4gICAgICAgICAgICBpZiAodGV4dE1lc2ggJiYgKGNpdHlTZXQuaW5kZXhPZihwdCkgIT0gLTEpKSB7XG4gICAgICAgICAgICAgICAgdmFyIHRleHRFZGdlID0gRWRnZS5jcmVhdGUocHQsIHBvaW50c1twdC5pZCArIDFdKTtcbiAgICAgICAgICAgICAgICBlZGdlcy5wdXNoKHRleHRFZGdlKTtcbiAgICAgICAgICAgICAgICBwdC5uZXh0cy5wdXNoKHRleHRFZGdlKTtcbiAgICAgICAgICAgIH1cblxuICAgICAgICB9XG4gICAgfSk7XG5cbiAgICByZXR1cm4gZWRnZXM7XG59OyIsIid1c2Ugc3RyaWN0JztcblxudmFyIHNxcnQgPSBNYXRoLnNxcnQ7XG52YXIgcG93ID0gTWF0aC5wb3c7XG52YXIgYWJzID0gTWF0aC5hYnM7XG52YXIgYXRhbiA9IE1hdGguYXRhbjtcblxudmFyIFZlY3RvciA9IHJlcXVpcmUoJy4vdmVjdG9yLmpzJyk7XG5cblxuZnVuY3Rpb24gRWRnZShwdEEsIHB0Qikge1xuICAgIHZhciBkaXN0YW5jZSA9IHNxcnQoIHBvdyhwdEEueCAtIHB0Qi54LCAyKSArIHBvdyhwdEEueSAtIHB0Qi55LCAyKSApO1xuXG4gICAgLy8gZmluZCBsaW5lIGVxdWF0aW9uIGF4ICsgYnkgKyBjID0gMFxuICAgIHZhciBhID0gMTtcbiAgICB2YXIgYiA9IC0gKHB0Qi54IC0gcHRBLngpIC8gKHB0Qi55IC0gcHRBLnkpO1xuXG4gICAgLy8gb3JpZW50YXRlIHZlY3RvciAoYSxiKVxuICAgIGlmIChiIDwgMCl7XG4gICAgICAgIGIgPSAtYjtcbiAgICAgICAgYSA9IC1hO1xuICAgIH1cblxuICAgIC8vIG5vcm1hbGl6ZSB2ZWN0b3IgKGEsYilcbiAgICB2YXIgbiA9IG5ldyBWZWN0b3IoYSwgYik7XG4gICAgbi5ub3JtYWxpemUoKTtcblxuICAgIHZhciBjID0gLSAoYSAqIHB0QS54ICsgYiAqIHB0QS55KTtcblxuICAgIC8vIC8vIGNhbGN1bGF0ZSB2ZWN0b3IgZGlyZWN0b3JcbiAgICB2YXIgdiA9IG5ldyBWZWN0b3IocHRCLnggLSBwdEEueCwgcHRCLnkgLSBwdEEueSk7XG4gICAgXG4gICAgdi5ub3JtYWxpemUoKTtcblxuXG4gICAgdGhpcy5pZCA9IHVuZGVmaW5lZDtcbiAgICB0aGlzLnB0MSA9IHB0QTtcbiAgICB0aGlzLnB0MiA9IHB0QjtcbiAgICB0aGlzLmRpcmVjdGlvbiA9IHY7XG4gICAgdGhpcy5vcnRoRGlyZWN0aW9uID0gbjsgXG4gICAgdGhpcy5kaXN0YW5jZSA9IGRpc3RhbmNlO1xuICAgIHRoaXMucGhlcm9tb24gPSAxL2Rpc3RhbmNlO1xuICAgIHRoaXMubGluZSA9IHtcbiAgICAgICAgYTogYSxcbiAgICAgICAgYjogYixcbiAgICAgICAgYzogYyxcbiAgICB9XG59XG5cblxuLy8gc3RhdGljIG1ldGhvZHNcbkVkZ2UuY3JlYXRlID0gZnVuY3Rpb24ocHRBLCBwdEIpIHtcbiAgICB2YXIgZWRnZSA9IG5ldyBFZGdlKHB0QSwgcHRCKTtcbiAgICByZXR1cm4gZWRnZTtcbn1cblxuXG4vLyBtZXRob2RzXG5FZGdlLnByb3RvdHlwZSA9IHtcblxuICAgIGdldE90aGVyUG9pbnQ6IGZ1bmN0aW9uKHBvaW50KSB7XG4gICAgICAgIGlmIChwb2ludCA9PSB0aGlzLnB0MSlcbiAgICAgICAgICAgIHJldHVybiB0aGlzLnB0MjtcbiAgICAgICAgZWxzZSBpZiAocG9pbnQgPT0gdGhpcy5wdDIpXG4gICAgICAgICAgICByZXR1cm4gdGhpcy5wdDE7XG4gICAgICAgIGVsc2VcbiAgICAgICAgICAgIGNvbnNvbGUubG9nKFwiRXJyb3JcIik7XG4gICAgfSxcblxuICAgIGNhbGN1bGF0ZURpc3RhbmNlOiBmdW5jdGlvbih4LCB5KSB7XG4gICAgICAgIHZhciBhID0gdGhpcy5saW5lLmEsXG4gICAgICAgICAgICBiID0gdGhpcy5saW5lLmIsXG4gICAgICAgICAgICBjID0gdGhpcy5saW5lLmM7XG4gICAgICAgIHJldHVybiBhYnMoYSAqIHggKyBiICogeSArIGMpIC8gTWF0aC5zcXJ0KE1hdGgucG93KGEsMikgKyBNYXRoLnBvdyhiLDIpKTtcbiAgICB9LFxuXG59XG5tb2R1bGUuZXhwb3J0cyA9IEVkZ2U7IiwiJ3VzZSBzdHJpY3QnXG5cbnZhciBhbnRGdW5jdGlvbiA9IHJlcXVpcmUoJy4vYW50LmpzJyk7XG5cbi8vIHZhciBOQkFOVFMgPSA0MDAwO1xuXG5tb2R1bGUuZXhwb3J0cyA9IGZ1bmN0aW9uIChjb250YWluZXIsIHBvaW50c0luZm9zLCBvcHRpb25zKSB7XG5cblx0dmFyIEFudCA9IGFudEZ1bmN0aW9uKGNvbnRhaW5lciwgcG9pbnRzSW5mb3MsIG9wdGlvbnMpO1xuXG5cdHZhciBwb3B1bGF0aW9uID0gbmV3IEFycmF5KG9wdGlvbnMubmJBbnRzKTtcblx0dmFyIHBvc3NpYmxlU3RhcnRQb2ludHNJZCA9IHBvaW50c0luZm9zLnBvc3NpYmxlU3RhcnRQb2ludHNJZDtcblxuXHRmb3IgKHZhciBpID0gMDsgaSA8IG9wdGlvbnMubmJBbnRzOyBpKyspIHtcblx0XHR2YXIgbmV3QW50ID0gbmV3IEFudChBbnQuZ2VuZXJhdGVSYW5kU3RhcnRQb2ludCgpKTtcblx0XHRuZXdBbnQuc2V0RGlyZWN0aW9uKCk7XG5cdFx0cG9wdWxhdGlvbltpXSA9IG5ld0FudDtcblx0fVxuXG5cdHJldHVybiBwb3B1bGF0aW9uO1xuXG59IiwiJ3VzZSBzdHJpY3QnXG5cbnZhciBwYXJzZSA9IHJlcXVpcmUoJ3BhcnNlLXN2Zy1wYXRoJyk7XG5cbnZhciByYW5nZSA9IHJlcXVpcmUoJy4vdXRpbGl0aWVzLmpzJykucmFuZ2U7XG5cbnZhciBQb2ludCA9IHJlcXVpcmUoJy4vcG9pbnQuanMnKTtcblxudmFyIHJhbmRvbSA9IE1hdGgucmFuZG9tO1xuXG4vLyB2YXIgbmJSYW5kb21Qb2ludHMgPSA1MDA7XG4vLyB2YXIgbmJTdGFydFBvaW50cyA9IDIwO1xuXG52YXIgbmJDaXR5ID0gMjtcblxudmFyIHRleHRNZXNoID0gdHJ1ZTtcblxuLy8gRnJhbWUgZGVmaW5pdGlvblxudmFyIHhJbml0ID0gMCwgeUluaXQgPSAwO1xudmFyIHcgPSAxLFxuICAgIGggPSAxO1xuXG52YXIgQWNoYXIgPSBcImMgNC42MDExLC0xMS43MTA0NyA5LjIwODM1LC0yMy40MjAwNiAxMy44MTk5LC0zNS4xMjg5OCA0LjYxMTU2LC0xMS43MDg5MiA5LjIyNzQxLC0yMy40MTcxOCAxMy44NDU3MywtMzUuMTI1IDQuNjE4MzEsLTExLjcwNzgyIDkuMjM5MDgsLTIzLjQxNTE5IDEzLjg2MDQ2LC0zNS4xMjIzMyA0LjYyMTM4LC0xMS43MDcxNCA5LjI0MzM2LC0yMy40MTQwNiAxMy44NjQxMSwtMzUuMTIwOTcgNC42MjA3NCwtMTEuNzA2OTEgOS4yNDAyNSwtMjMuNDEzODEgMTMuODU2NjcsLTM1LjEyMDkyIDQuNjE2NDEsLTExLjcwNzEyIDkuMjI5NzQsLTIzLjQxNDQ0IDEzLjgzODE0LC0zNS4xMjIxOCA0LjYwODM5LC0xMS43MDc3NSA5LjIxMTg2LC0yMy40MTU5MiAxMy44MDg1NCwtMzUuMTI0NzQgNC41OTY2OCwtMTEuNzA4ODEgOS4xODY1OCwtMjMuNDE4Mjc3IDEzLjc2Nzg1LC0zNS4xMjg2MDMgMTIuOTk5MjMsLTMuMzU3ODU1IDI0LjMwMDY5LC02LjgyMTE0NCAzMy44Njg5MywtNC40NjQxMSA5LjU2ODI1LDIuMzU3MDM1IDE3LjQwMzI4LDEwLjUzNDM5MyAyMy40Njk2NywzMC40NTc4MzMgNC4xNzU5OCwxMC42MjExNCA4LjM2MDMxLDIxLjIzODgyIDEyLjU0OTQyLDMxLjg1NDUzIDQuMTg5MSwxMC42MTU3MSA4LjM4Mjk4LDIxLjIyOTQzIDEyLjU3ODA3LDMxLjg0MjY2IDQuMTk1MDksMTAuNjEzMjMgOC4zOTE0LDIxLjIyNTk1IDEyLjU4NTM1LDMxLjgzOTY1IDQuMTkzOTUsMTAuNjEzNjkgOC4zODU1NSwyMS4yMjgzNiAxMi41NzEyMywzMS44NDU0OSA0LjE4NTY4LDEwLjYxNzEyIDguMzY1NDQsMjEuMjM2NyAxMi41MzU3MiwzMS44NjAyMSA0LjE3MDI4LDEwLjYyMzUgOC4zMzEwOCwyMS4yNTA5MyAxMi40Nzg4MywzMS44ODM3NyA0LjE0Nzc1LDEwLjYzMjgzIDguMjgyNDUsMjEuMjcxMDcgMTIuNDAwNTMsMzEuOTE2MiA0LjExODA4LDEwLjY0NTEyIDguMjE5NTUsMjEuMjk3MTIgMTIuMzAwODUsMzEuOTU3NDkgLTUuOTA4OTYsNi45NTU2MSAtMjQuNjE2MTcsMS4xMTI5OCAtMzUuNTkzNzIsMyAtMTYuNzUyNDgsMi44NDg1OSAtMjYuOTY0MjEsLTAuNDE0MTYgLTI4LjQwNjI4LC0xOSAtMy4xNzcyNiwtOC40MjE1NyAtNi4zNTYwNiwtMTYuODc5MjQgLTkuNDc5OTcsLTI1LjM0ODU0IC0zLjEyMzkxLC04LjQ2OTI5IC02LjE5MjkzLC0xNi45NTAyIC05LjE1MDYyLC0yNS40MTgyNiAtMTYuNDY3MjMsLTAuODAwNTMgLTM1LjE3NzY4LC0yLjc0MjgxIC01My4zMzcyNywtMy4xMTQzOSAtMTguMTU5NiwtMC4zNzE1OCAtMzUuNzY4MzQsMC44Mjc1MyAtNTAuMDMyMTQsNi4zMDk3NiAtNC4xNTY3OSwxMC45MzEyNCAtOC4xMzgwNiwyMS45MjY5IC0xMi4xNTc4NSwzMi45MDU0NiAtNC4wMTk4LDEwLjk3ODU1IC04LjA3ODEzLDIxLjk0MDAxIC0xMi4zODkwMywzMi44MDI4MiAtOS41Mjc1OCwwLjgyNzQ3IC0xOS4xMDQ1NywwLjk2NDIzIC0yOC42OTI4MiwwLjkzMzYzIC05LjU4ODI0LC0wLjAzMDYgLTE5LjE4NzczLC0wLjIyODU1IC0yOC43NjAzLC0wLjA3MDUgbCAwLC0yLjE5NzkxIHogbSAxNzQsLTEwOSBjIC0zLjE3NDYyLC05LjUxMTM2IC02LjQ0ODM1LC0xOC45OTE3NCAtOS43NTQzLC0yOC40NjIyNCAtNi44MTU1OSwtMTkuNTE0MTIgLTYuMjQyMDIsLTI1LjM3OTQ2IC0xMi4zMTc0LC00My42MDc4OCAtMy4xNjcyMiwtOS41MTU0MyAtMTMuNjAyNTcsLTMyLjMyODY2IC0xNi40OTk3MywtNDEuOTI5ODggLTYuMDg5ODksMTIuODg1NzYgLTExLjEzMiwyNi41OTI3MiAtMTUuOTM0MTUsNDAuNDc0NzYgLTQuODAyMTUsMTMuODgyMDMgLTkuMzY0MzUsMjcuOTM5MTUgLTE0LjQ5NDQyLDQxLjUyNTI0IC0yLjM5NDI1LDEyLjA1ODAxIC0yMi44ODE1LDQwLjExMTE2IDIuMTg3NDUsMzQgMTEuMDUyNTYsLTAuNDQ0ODYgMjIuNTI5MzksMC4xMTA2MSAzMy44NTYyNCwwLjI0OTU1IDExLjMyNjg0LDAuMTM4OTUgMjIuNTAzNjksLTAuMTM4NjMgMzIuOTU2MzEsLTIuMjQ5NTUgeiBcIjtcbnZhciBOY2hhciA9IFwiYyAwLDAgMCwtMjMuODMzMzQgMCwtMzUuNzUgMCwtMTEuOTE2NjcgMCwtMjMuODMzMzQgMCwtMzUuNzUgMCwtMTEuOTE2NjYgMCwtMjMuODMzMzMgMCwtMzUuNzUgMCwtMTEuOTE2NjY3IDAsLTIzLjgzMzMzNSAwLC0zNS43NTAwMDMgMTQuNjQ3MjksMS4wNTMxMyAzMS4wMzE5MywtMi4wODgwOCA0NC42MDY3OSwxLjU1NDE1IDYuMTgwNzYsNy42MTA5OTYgMTIuMzU0MzgsMTUuMjMxMzQyIDE4LjUxOTczLDIyLjg2MTAxMyA2LjE2NTM1LDcuNjI5NjcgMTIuMzIyNDQsMTUuMjY4NjYgMTguNDcwMTQsMjIuOTE2OTYgNi4xNDc3LDcuNjQ4MjkgMTIuMjg2MDIsMTUuMzA1ODggMTguNDEzODMsMjIuOTcyNzQgNi4xMjc4MSw3LjY2Njg2IDEyLjI0NTEyLDE1LjM0MyAxOC4zNTA4MSwyMy4wMjgzOCA2LjEwNTY4LDcuNjg1MzggMTIuMTk5NzUsMTUuMzgwMDEgMTguMjgxMDcsMjMuMDgzODYgNi4wODEzMiw3LjcwMzg0IDEyLjE0OTksMTUuNDE2OTEgMTguMjA0NjIsMjMuMTM5MTggNi4wNTQ3Miw3LjcyMjI2IDEyLjA5NTU3LDE1LjQ1MzcyIDE4LjEyMTQ1LDIzLjE5NDM1IDYuMDI1ODcsNy43NDA2MyAxMi4wMzY3NiwxNS40OTA0MyAxOC4wMzE1NiwyMy4yNDkzNyAwLjQ1NTgsLTE1LjQ5MTQgMC43MjY1NSwtMzAuOTg3MjUgMC44ODExOSwtNDYuNDg1ODIgMC4xNTQ2NCwtMTUuNDk4NTggMC4xOTMxOSwtMzAuOTk5ODggMC4xODQ1OCwtNDYuNTAyMTggLTAuMDA5LC0xNS41MDIzIC0wLjA2NDMsLTMxLjAwNTYgLTAuMDk4MywtNDYuNTA4MTcgLTAuMDM0LC0xNS41MDI1OCAtMC4wNDYxLC0zMS4wMDQ0MzIgMC4wMzI1LC00Ni41MDM4MzMgOS4zMzMzMywwIDE4LjY2NjY3LDAgMjgsMCA5LjMzMzMzLDAgMTguNjY2NjYsMCAyOCwwIDAsMTEuOTE2NjY3IDAsMjMuODMzMzMyIDAsMzUuNzUwMDAzIDAsMTEuOTE2NjYgMCwyMy44MzMzMyAwLDM1Ljc1IDAsMTEuOTE2NjYgMCwyMy44MzMzMyAwLDM1Ljc1IDAsMTEuOTE2NjcgMCwyMy44MzMzMyAwLDM1Ljc1IDAsMTEuOTE2NjYgMCwyMy44MzMzMyAwLDM1Ljc1IDAsMTEuOTE2NjcgMCwyMy44MzMzNCAwLDM1Ljc1IDAsMTEuOTE2NjYgMCwyMy44MzMzMyAwLDM1Ljc1IDAsMTEuOTE2NjcgMCwyMy44MzMzNCAwLDM1Ljc1IC0xNC42NDQxLC0xLjA1MzQ4IC0zMS4wMjYsMi4wODg0NyAtNDQuNTk3NjQsLTEuNTU0MTYgLTYuMjcyNDcsLTcuMzgwMzcgLTEyLjQ4NDM2LC0xNC44MTY5MiAtMTguNjUxMjMsLTIyLjI5NTA3IC02LjE2Njg4LC03LjQ3ODE1IC0xMi4yODg3MywtMTQuOTk3OSAtMTguMzgxMTEsLTIyLjU0NDY1IC02LjA5MjM4LC03LjU0Njc1IC0xMi4xNTUyOCwtMTUuMTIwNTEgLTE4LjIwNDI1LC0yMi43MDY2NyAtNi4wNDg5NiwtNy41ODYxNyAtMTIuMDgzOTksLTE1LjE4NDc2IC0xOC4xMjA2MywtMjIuNzgxMTYgLTYuMDM2NjQsLTcuNTk2NCAtMTIuMDc0ODksLTE1LjE5MDYyIC0xOC4xMzAzLC0yMi43NjgwNyAtNi4wNTU0MSwtNy41Nzc0NSAtMTIuMTI3OTcsLTE1LjEzODEzIC0xOC4yMzMyMywtMjIuNjY3NDQgLTYuMTA1MjYsLTcuNTI5MzIgLTEyLjI0MzIyLC0xNS4wMjcyNyAtMTguNDI5NDIsLTIyLjQ3OTI2IC02LjE4NjE5LC03LjQ1MTk5IC0xMi40MjA2MywtMTQuODU4MDMgLTE4LjcxODg2LC0yMi4yMDM1MiAtMC4yMjc5MSwxNS4xNjQ5NiAtMC4zNjY3NCwzMC4zMzA4IC0wLjQ0ODksNDUuNDk3MTggLTAuMDgyMiwxNS4xNjYzNyAtMC4xMDc2NiwzMC4zMzMyOSAtMC4xMDg4OSw0NS41MDA0IC0wLjAwMSwxNS4xNjcxMiAwLjAyMTgsMzAuMzM0NDMgMC4wMzY3LDQ1LjUwMTYyIDAuMDE0OSwxNS4xNjcxOCAwLjAyMTYsMzAuMzM0MjIgLTAuMDEyMiw0NS41MDA4IC05LjMzMzMzLDAgLTE4LjY2NjY3LDAgLTI4LDAgLTkuMzMzMzMsMCAtMTguNjY2NjYsMCAtMjgsMCAwLC0xMS45MTY2NyAwLC0yMy44MzMzNCAwLC0zNS43NSAwLC0xMS45MTY2NyAwLC0yMy44MzMzNCAwLC0zNS43NSAwLC0xMS45MTY2NiAwLC0yMy44MzMzMyAwLC0zNS43NSAwLC0xMS45MTY2NyAwLC0zNS43NSAtMTBlLTYsLTM1Ljc1IHogXCI7XG52YXIgVGNoYXIgPSBcImMgMCwtOS45MTY2NyAwLC0xOS44MzMzNCAwLC0yOS43NSAwLC05LjkxNjY3IDAsLTE5LjgzMzM0IDAsLTI5Ljc1IDAsLTkuOTE2NjYgMCwtMTkuODMzMzMgMCwtMjkuNzUgMCwtOS45MTY2NiAtMS40NzM5NCwtMTkuODMzMzMgLTEuNDczOTQsLTI5Ljc1IDAsLTE1LjMzMzM0IC0yOS4xOTI3MywwIC00NC41MjYwNiwwIC0xNS4zMzMzMywwIC0zMC42NjY2NywwIC00NiwwIDAsLTggMCwtMTYgMCwtMjQuMDAwMDAyIDAsLTggMCwtMTYuMDAwMDAxIDAsLTI0LjAwMDAwMSA5LjkxNjY2LDAgMTkuODMzMzMsMCAyOS43NSwwIDkuOTE2NjcsMCAxOS44MzMzMywwIDI5Ljc1LDAgOS45MTY2NiwwIDE5LjgzMzMzLDAgMjkuNzUsMCA5LjkxNjY3LDAgMTkuODMzMzMsMCAyOS43NSwwIDkuOTE2NjYsMCAxOS44MzMzMywwIDI5Ljc1LDAgOS45MTY2NywwIDE5LjgzMzMzLDAgMjkuNzUsMCA5LjkxNjY2LDAgMTkuODMzMzMsMCAyOS43NSwwIDkuOTE2NjcsMCAxOS44MzMzMywwIDI5Ljc1LDAgMCw4IDAsMTYuMDAwMDAxIDAsMjQuMDAwMDAxIDAsOC4wMDAwMDIgMCwxNi4wMDAwMDIgMCwyNC4wMDAwMDIgLTE1LDAgLTMwLDAgLTQ1LDAgLTE1LDAgLTQ0Ljg4ODA5LC0xMy41MjU2NSAtNDUsMS40NzM5NCAtMC4wNzM1LDkuODU1ODggLTAuMTA2MTMsMTguMjM5OSAtMC4xMTI0OSwyOC4wOTkyMyAtMC4wMDYsOS44NTkzMiAwLjAxMzUsMTkuNzIwMDEgMC4wNDUsMjkuNTgxMzMgMC4wMzE1LDkuODYxMzIgMC4wNzQ1LDE5LjcyMzI5IDAuMTE0MywyOS41ODUxNyAwLjAzOTksOS44NjE4OCAwLjA3NjUsMTkuNzIzNjcgMC4wOTU0LDI5LjU4NDY3IDAuMDE4OSw5Ljg2MDk5IDAuMDE5OSwxOS43MjExOCAtMC4wMTE3LDI5LjU3OTg0IC0wLjAzMTUsOS44NTg2NiAtMC4wOTU2LDE5LjcxNTggLTAuMjA2ODgsMjkuNTcwNjkgLTAuMTExMzEsOS44NTQ4OCAtMC4yNjk4NSwxOS43MDc1MiAtMC40OTAzMywyOS41NTcxOSAtMC4yMjA0Nyw5Ljg0OTY3IC0wLjUwMjg5LDE5LjY5NjM2IC0wLjg2MTkzLDI5LjUzOTM3IC03LjM5NjI0LC0wLjk5OTY0IC0xOC45NTY1OCwwLjk2NzI2IC0yOS43MDkxMiwxLjgyOTcxIC0xMC43NTI1MywwLjg2MjQ1IC0yNC40ODU1NiwyLjAyNjA5IC0yNC44NjIzMSwtNC43OTY5NiAtMC42NjIyMywtMTEuOTkzMTQgLTAuNjYyMjMsLTIxLjU0MzQ5IC0wLjQ5NjY3LC0zMC40ODMxNCAwLjE2NTU2LC04LjkzOTY1IDAuNDk2NjcsLTE3LjI2ODYgMC40OTY2NywtMjYuODE4OTUgMCwtOS41NTAzNSAwLC0xOS4xMDA3IDAsLTI4LjY1MTA1IDAsLTkuNTUwMzUgMCwtMTkuMTAwNjkgMCwtMjguNjUxMDQgelwiO1xudmFyIFNjaGFyID0gXCJjIC0xNi42NTg3LC0yLjM0Mzg3IC0zMy4zNTk5NSwtNi4yMzUxNSAtNDkuMjUxMzcsLTEyLjA4ODMxIC0xNS44OTE0MiwtNS44NTMxNiAtMzAuOTczMDUsLTEzLjY2ODIgLTQ0LjM5MjYsLTIzLjg1OTU3IDMuNzM3MTUsLTguMDI2MzkgNy45NjA2NSwtMTUuNzkzNzEgMTIuMTkxMDEsLTIzLjU1MjI4IDQuMjMwMzcsLTcuNzU4NTggOC40Njc1OSwtMTUuNTA4NDIgMTIuMjMyMTksLTIzLjQ5OTg0IDguODI5MDIsNi42OTE2MiAxOC40ODIzNSwxMi43NTY2MyAyOC42NzQxMiwxNy44NjM0NCAxMC4xOTE3OCw1LjEwNjgxIDIwLjkyMjA1LDkuMjU1NDIgMzEuOTA0ODUsMTIuMTE0MjEgMTAuOTgyOSwyLjg1ODggMjIuMjE4Myw0LjQyNzggMzMuNDIwNiw0LjM3NTM5IDExLjIwMjIsLTAuMDUyNCAyMi4zNzEzLC0xLjcyNjIyIDMzLjIyMTIsLTUuMzUzMDQgMTQuMjgwNCwtNi45MjE2OSAxNy4zMDMzLC0xOS42NDQzNiAxMy40OTQ2LC0zMS4xNTAyOCAtMy44MDg3LC0xMS41MDU5MSAtMTQuNDQ4OSwtMjEuNzk1MDcgLTI3LjQ5NDYsLTIzLjg0OTcyIC0xMC4xNjIsLTQuMDc2ODYgLTIxLjEwNTIsLTcuMTQ4NjMgLTMyLjE5NzUsLTEwLjExMDk3IC0xMS4wOTI0LC0yLjk2MjM1IC0yMi4zMzQsLTUuODE1MjYgLTMzLjA5MzEsLTkuNDU0MzkgLTEwLjc1OTEsLTMuNjM5MTQgLTIxLjAzNTYyLC04LjA2NDQ5IC0zMC4xOTc3OSwtMTQuMTcxNzEgLTkuMTYyMTYsLTYuMTA3MjMgLTE3LjIwOTk2LC0xMy44OTYzMiAtMjMuNTExNTgsLTI0LjI2MjkzIC01LjE3OTQsLTEwLjk3MDU4IC03LjQ4NTQ0LC0yMi43NzU4MyAtNy4zMDUxMiwtMzQuNTI5NjYgMC4xODAzMSwtMTEuNzUzODIgMi44NDY5OCwtMjMuNDU2MiA3LjYxMjk5LC0zNC4yMjEwMyA0Ljc2NjAxLC0xMC43NjQ4MyAxMS42MzEzNywtMjAuNTkyMSAyMC4yMDkwNSwtMjguNTk1NjkyIDguNTc3NjksLTguMDAzNTkyIDE4Ljg2NzcxLC0xNC4xODM1MDYgMzAuNDgzMDUsLTE3LjY1MzYyMSAxMi4xODIxLC00LjAyNDE0NSAyNC44Nzc3LC02LjM1NzAwOSAzNy42OTg4LC03LjEyMTkwMSAxMi44MjExLC0wLjc2NDg5MiAyNS43Njc3LDAuMDM4MTkgMzguNDUxNywyLjI4NTkzMyAxMi42ODQxLDIuMjQ3NzQ0IDI1LjEwNTUsNS45NDAxNTMgMzYuODc2NSwxMC45NTM5MTggMTEuNzcwOSw1LjAxMzc2NCAyMi44OTEyLDExLjM0ODg4NSAzMi45NzMsMTguODgyMDUzIC00LjQ3NDUsNy4wMDM2OCAtOC4yNTgxLDE0LjcwNzA4IC0xMi4xODg2LDIyLjIzNDk5IC0zLjkzMDYsNy41Mjc5IC04LjAwOCwxNC44ODAzMSAtMTMuMDcwMSwyMS4xODE5OCAtMTUuMDcyOCwtMTEuMTUxNTggLTMzLjUwMDksLTIwLjU0NDkxIC01Mi42MTE0LC0yNC45MTcyNiAtMTkuMTEwNSwtNC4zNzIzNSAtMzguOTAzNCwtMy43MjM3MiAtNTYuNzA1OCw1LjIwODYxIC0xMC4wMzIxLDcuMDI3ODkgLTEzLjMzNDgsMTguODQ2OTUgLTExLjE1NjEsMjkuNTM1OTYgMi4xNzg3LDEwLjY4OTAyIDkuODM4NywyMC4yNDc5OSAyMS43MzIsMjIuNzU1NzIgMTAuMTMzNSw0LjQzMzA3IDIxLjAyNDEsNy42OTQwNyAzMi4xMTAxLDEwLjY5MDIzIDExLjA4NjEsMi45OTYxNiAyMi4zNjc2LDUuNzI3NDcgMzMuMjgzLDkuMTAxMTcgMTAuOTE1NSwzLjM3MzcgMjEuNDY0OCw3LjM4OTc4IDMxLjA4NjUsMTIuOTU1NDcgOS42MjE3LDUuNTY1NjkgMTguMzE1NywxMi42ODA5OSAyNS41MjA0LDIyLjI1MzEzIDYuNjMwMSw5LjcxNjM1IDEwLjU1ODcsMjAuODgwMjMgMTIuMDA2NCwzMi40MTIzMSAxLjQ0NzYsMTEuNTMyMDcgMC40MTQyLDIzLjQzMjM0IC0yLjg3OTcsMzQuNjIxNDYgLTMuMjkzOSwxMS4xODkxMiAtOC44NDg0LDIxLjY2NzA5IC0xNi40NDMsMzAuMzU0NTcgLTcuNTk0Niw4LjY4NzQ5IC0xNy4yMjkzLDE1LjU4NDQ5IC0yOC42ODM3LDE5LjYxMTY2IC0xMy4wMDkxLDUuOTI5MjggLTI3LjAxOTcsOC42Mzg2OSAtNDEuMjcyOCw5LjYzNjA4IC0xNC4yNTMsMC45OTczOCAtMjguNzQ4NCwwLjI4Mjc0IC00Mi43MjcyLC0wLjYzNjA4IHogXCI7XG5cbnZhciBzdmdTdHJpbmcgPSBcIm0gMCwwIFwiICsgQWNoYXIgK1wibSAxMjYsLTMxIFwiICsgTmNoYXIgKyBcIm0gMzc2LDI0IFwiKyBUY2hhciArXCJtIDI1MCwxMjAgXCIrIFNjaGFyO1xuXG4gXG5mdW5jdGlvbiBzdmdUb1BvaW50cyhzdmdTdHJpbmcpIHtcbiAgICB2YXIgcG9pbnRzID0gW107XG4gICAgdmFyIGVkZ2VzID0gT2JqZWN0LmNyZWF0ZShudWxsKTtcblxuICAgIHZhciBiZWdpbmluZ1BhdGg7XG5cbiAgICB2YXIgWCA9IDA7XG4gICAgdmFyIFkgPSAwO1xuICAgIHZhciBuYlBvaW50cyA9IDA7XG4gICAgdmFyIHByZXZQb2ludDtcblxuICAgIHZhciBjb21tYW5kcyA9IHBhcnNlKHN2Z1N0cmluZylcbiAgICBmb3IgKHZhciBpPTA7IGk8Y29tbWFuZHMubGVuZ3RoOyBpKyspe1xuICAgICAgICB2YXIgY29tbWFuZCA9IGNvbW1hbmRzW2ldO1xuICAgICAgICBzd2l0Y2ggKGNvbW1hbmRbMF0pIHtcbiAgICAgICAgICAgIGNhc2UgXCJtXCI6XG4gICAgICAgICAgICAgICAgWCArPSBjb21tYW5kWzFdO1xuICAgICAgICAgICAgICAgIFkgKz0gY29tbWFuZFsyXTtcbiAgICAgICAgICAgICAgICBwcmV2UG9pbnQgPSB1bmRlZmluZWQ7XG4gICAgICAgICAgICAgICAgYmVnaW5pbmdQYXRoID0gbmJQb2ludHM7XG4gICAgICAgICAgICAgICAgYnJlYWs7XG4gICAgICAgICAgICBjYXNlIFwiTVwiOlxuICAgICAgICAgICAgICAgIFggPSBjb21tYW5kWzFdO1xuICAgICAgICAgICAgICAgIFkgPSBjb21tYW5kWzJdO1xuICAgICAgICAgICAgICAgIHByZXZQb2ludCA9IHVuZGVmaW5lZDtcbiAgICAgICAgICAgICAgICBiZWdpbmluZ1BhdGggPSBuYlBvaW50cztcbiAgICAgICAgICAgICAgICBicmVhazsgIFxuICAgICAgICAgICAgY2FzZSBcImNcIjpcbiAgICAgICAgICAgICAgICBYICs9IGNvbW1hbmRbNV07XG4gICAgICAgICAgICAgICAgWSArPSBjb21tYW5kWzZdO1xuICAgICAgICAgICAgICAgIHBvaW50cy5wdXNoKHtpZDpuYlBvaW50cywgeDpYLCB5Oll9KTtcblxuICAgICAgICAgICAgICAgIGlmIChwcmV2UG9pbnQgIT0gdW5kZWZpbmVkKSB7XG4gICAgICAgICAgICAgICAgICAgIGVkZ2VzW3ByZXZQb2ludF0gPSBuYlBvaW50cztcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgcHJldlBvaW50ID0gbmJQb2ludHM7XG4gICAgICAgICAgICAgICAgbmJQb2ludHMrKztcbiAgICAgICAgICAgICAgICBicmVhaztcbiAgICAgICAgICAgIGNhc2UgXCJ6XCI6XG4gICAgICAgICAgICAgICAgZWRnZXNbcHJldlBvaW50XSA9IG5iUG9pbnRzO1xuICAgICAgICAgICAgICAgIGJlZ2luaW5nUGF0aCA9IHVuZGVmaW5lZDtcbiAgICAgICAgICAgICAgICBwcmV2UG9pbnQgPSB1bmRlZmluZWQ7XG4gICAgICAgICAgICAgICAgYnJlYWs7ICAgIFxuICAgICAgICB9XG4gICAgfVxuICAgIHJldHVybiB7cG9pbnRzIDogcG9pbnRzLCBlZGdlcyA6IGVkZ2VzfTtcbn1cblxuLy8gaW5pdGlhbGl6ZSBwb2ludHNcblxubW9kdWxlLmV4cG9ydHMgPSBmdW5jdGlvbihuYlN0YXJ0UG9pbnRzLCBuYlJhbmRvbVBvaW50cyl7XG4gICAgdmFyIHBvaW50cyA9IFtdO1xuICAgIHZhciBmb3JjZWRFZGdlcztcbiAgICB2YXIgY2l0eVNldDtcblxuICAgIGlmICh0ZXh0TWVzaCl7XG5cbiAgICAgICAgdmFyIG15VGV4dCA9IHN2Z1RvUG9pbnRzKHN2Z1N0cmluZyk7XG4gICAgICAgIHBvaW50cyA9IG15VGV4dC5wb2ludHM7XG4gICAgICAgIGZvcmNlZEVkZ2VzID0gbXlUZXh0LmVkZ2VzO1xuICAgICAgICBjaXR5U2V0ID0gcmFuZ2UoMCwgcG9pbnRzLmxlbmd0aCk7XG5cbiAgICAgICAgdmFyIHNjYWxlWCA9IDAuNTtcbiAgICAgICAgdmFyIHNjYWxlWSA9IDAuMjtcbiAgICAgICAgdmFyIGRlbHRhWCA9IDAuMjU7XG4gICAgICAgIHZhciBkZWx0YVkgPSAwLjM7XG5cbiAgICAgICAgLy8gc2NhbGUgcG9pbnRzIHRvIFswLDFdICsgZGVsdGFcbiAgICAgICAgdmFyIG1heFggPSBNYXRoLm1heC5hcHBseShNYXRoLCBwb2ludHMubWFwKGZ1bmN0aW9uKHApe3JldHVybiBwLnh9KSk7XG4gICAgICAgIHZhciBtaW5YID0gTWF0aC5taW4uYXBwbHkoTWF0aCwgcG9pbnRzLm1hcChmdW5jdGlvbihwKXtyZXR1cm4gcC54fSkpO1xuICAgICAgICB2YXIgbWF4WSA9IE1hdGgubWF4LmFwcGx5KE1hdGgsIHBvaW50cy5tYXAoZnVuY3Rpb24ocCl7cmV0dXJuIHAueX0pKTtcbiAgICAgICAgdmFyIG1pblkgPSBNYXRoLm1pbi5hcHBseShNYXRoLCBwb2ludHMubWFwKGZ1bmN0aW9uKHApe3JldHVybiBwLnl9KSk7XG4gICAgICAgIHBvaW50cyA9IHBvaW50cy5tYXAoZnVuY3Rpb24ocCl7XG4gICAgICAgICAgICB2YXIgeCA9IHNjYWxlWCAqIChwLngtbWluWCkvKG1heFgtbWluWCkgKyBkZWx0YVg7XG4gICAgICAgICAgICB2YXIgeSA9IHNjYWxlWSAqIChwLnktbWluWSkvKG1heFktbWluWSkgKyBkZWx0YVk7XG4gICAgICAgICAgICB2YXIgbmV3UG9pbnQgPSBuZXcgUG9pbnQoeCwgeSk7XG4gICAgICAgICAgICBuZXdQb2ludC5pZCA9IHAuaWQ7XG5cbiAgICAgICAgICAgIHJldHVybiBuZXdQb2ludDtcbiAgICAgICAgfSk7XG5cbiAgICAgICAgLy8gb25seSBhZGQgcmFuZG9tIHBvaW50c1xuICAgICAgICB2YXIgbmJQb2ludHMgPSBwb2ludHMubGVuZ3RoO1xuICAgICAgICBmb3IodmFyIGk9MDsgaTxuYlJhbmRvbVBvaW50czsgKytpKSB7XG5cbiAgICAgICAgICAgIHZhciB4ID0gcmFuZG9tKCk7XG4gICAgICAgICAgICB2YXIgeSA9IHJhbmRvbSgpO1xuXG4gICAgICAgICAgICB2YXIgbmV3UG9pbnQgPSBuZXcgUG9pbnQoeCwgeSk7XG4gICAgICAgICAgICBuZXdQb2ludC5pZCA9IG5iUG9pbnRzO1xuXG4gICAgICAgICAgICBwb2ludHMucHVzaChuZXdQb2ludCk7XG5cbiAgICAgICAgICAgIG5iUG9pbnRzKys7XG4gICAgICAgIH1cblxuICAgIH0gZWxzZSB7XG4gICAgICAgIC8vYWRkIHJhbmRvbSBwb2ludHNcblxuICAgICAgICB2YXIgbmJQb2ludHMgPSAwO1xuICAgICAgICBmb3IodmFyIGk9MDsgaTxuYlJhbmRvbVBvaW50czsgKytpKSB7XG5cbiAgICAgICAgICAgIHZhciB4ID0gcmFuZG9tKCk7XG4gICAgICAgICAgICB2YXIgeSA9IHJhbmRvbSgpO1xuXG4gICAgICAgICAgICB2YXIgbmV3UG9pbnQgPSBuZXcgUG9pbnQoeCwgeSk7XG4gICAgICAgICAgICBuZXdQb2ludC5pZCA9IG5iUG9pbnRzO1xuXG4gICAgICAgICAgICBwb2ludHMucHVzaChuZXdQb2ludCk7XG4gICAgICAgICAgICBcbiAgICAgICAgICAgIG5iUG9pbnRzKys7XG4gICAgICAgIH1cblxuICAgICAgICBjaXR5U2V0ID0gcmFuZ2UoMCwgbmJDaXR5KTtcbiAgICAgICAgY29uc29sZS5sb2coY2l0eVNldCk7XG4gICAgfVxuXG5cbiAgICAvLyBpbml0aWFsaXplIHN0YXJ0IHBvaW50c1xuICAgIHZhciBwb3NzaWJsZVN0YXJ0UG9pbnRzSWQgPSBbXTtcblxuICAgIGlmIChuYlJhbmRvbVBvaW50cyAhPT0gMCkge1xuICAgICAgICBmb3IgKHZhciBpID0gMDsgaSA8IG5iU3RhcnRQb2ludHM7IGkrKyl7XG4gICAgICAgICAgICBwb3NzaWJsZVN0YXJ0UG9pbnRzSWQucHVzaChNYXRoLmZsb29yKG5iUmFuZG9tUG9pbnRzICogcmFuZG9tKCkpKTtcbiAgICAgICAgfVxuICAgIH0gZWxzZSB7XG4gICAgICAgIGZvciAodmFyIGkgPSAwOyBpIDwgbmJTdGFydFBvaW50czsgaSsrKXtcbiAgICAgICAgICAgIHBvc3NpYmxlU3RhcnRQb2ludHNJZC5wdXNoKE1hdGguZmxvb3IobmJQb2ludHMgKiByYW5kb20oKSkpO1xuICAgICAgICB9XG4gICAgfVxuICAgIFxuXG4gICAgcmV0dXJuIHtcbiAgICAgICAgdGV4dE1lc2g6IHRleHRNZXNoLFxuICAgICAgICBwb2ludHM6IHBvaW50cyxcbiAgICAgICAgY2l0eVNldDogY2l0eVNldCxcbiAgICAgICAgcG9zc2libGVTdGFydFBvaW50c0lkOiBwb3NzaWJsZVN0YXJ0UG9pbnRzSWQsXG4gICAgICAgIG5iUmFuZG9tUG9pbnRzOiBuYlJhbmRvbVBvaW50cyxcbiAgICAgICAgZm9yY2VkRWRnZXM6IGZvcmNlZEVkZ2VzXG4gICAgfTtcbn1cbiIsIid1c2Ugc3RyaWN0J1xuXG5tb2R1bGUuZXhwb3J0cyA9IGZ1bmN0aW9uIChjb250YWluZXIpe1xuXG5cdHZhciBtb3VzZSA9IHtcblx0ICAgIHg6IDAsXG5cdCAgICB5OiAwXG5cdH07XG5cblx0Y29udGFpbmVyLmFkZEV2ZW50TGlzdGVuZXIoICdtb3VzZW1vdmUnLCBmdW5jdGlvbihlKXtcblx0ICAgIHZhciByZWN0ID0gY29udGFpbmVyLmdldEJvdW5kaW5nQ2xpZW50UmVjdCgpO1xuXG5cdCAgICBtb3VzZS54ID0gKGUuY2xpZW50WCAtIHJlY3QubGVmdCApIC8gcmVjdC53aWR0aDtcblx0ICAgIG1vdXNlLnkgPSAoZS5jbGllbnRZIC0gcmVjdC50b3AgKS8gcmVjdC5oZWlnaHQ7XG5cdH0pO1xuXG5cdHJldHVybiBtb3VzZTtcblxufTtcbiIsIid1c2Ugc3RyaWN0J1xuXG5mdW5jdGlvbiBQb2ludCh4LCB5KSB7XG4gICAgdGhpcy5pZCA9IHVuZGVmaW5lZDsgICAgICAgICAgICAgICAgXG4gICAgdGhpcy54ID0geDtcbiAgICB0aGlzLnkgPSB5O1xuICAgIHRoaXMubmV4dHMgPSBbXTtcbn1cblxubW9kdWxlLmV4cG9ydHMgPSBQb2ludDsiLCIndXNlIHN0cmljdCdcblxudmFyIGFudEZ1bmN0aW9uID0gcmVxdWlyZSgnLi9hbnQuanMnKTtcblxudmFyIHJhbmRvbSA9IE1hdGgucmFuZG9tO1xuXG52YXIgUkFORE9NTVZUID0gMC4wMDM7XG52YXIgQU5UU0laRSA9IDAuMDAyO1xuXG5tb2R1bGUuZXhwb3J0cyA9IGZ1bmN0aW9uKGNvbnRhaW5lciwgaW5pdFZhcil7XG5cdFxuXHRpZighY29udGFpbmVyKVxuXHRcdHRocm93IG5ldyBUeXBlRXJyb3IoJ01pc3NpbmcgY29udGFpbmVyJyk7XG5cblx0dmFyIGVkZ2VzID0gaW5pdFZhci5lZGdlcztcblx0dmFyIHBvcHVsYXRpb24gPSBpbml0VmFyLnBvcHVsYXRpb247XG5cdHZhciBwb2ludHNJbmZvcyA9IGluaXRWYXIucG9pbnRzSW5mb3M7XG5cdHZhciBuYkFudHMgPSBwb3B1bGF0aW9uLmxlbmd0aDtcblxuXHR2YXIgYW5pbUlEO1xuXG5cdHZhciBjYW52YXNMaXN0ID0gZG9jdW1lbnQuZ2V0RWxlbWVudHNCeVRhZ05hbWUoXCJjYW52YXNcIik7XG5cdFxuXHRpZiAoY2FudmFzTGlzdC5sZW5ndGggPT09IDApe1xuXHRcdHZhciBjYW52YXMgPSBkb2N1bWVudC5jcmVhdGVFbGVtZW50KFwiY2FudmFzXCIpO1xuXHRcdHZhciByZWN0ID0gY29udGFpbmVyLmdldEJvdW5kaW5nQ2xpZW50UmVjdCgpO1xuXHRcdGNhbnZhcy53aWR0aCA9IHJlY3Qud2lkdGg7XG5cdFx0Y2FudmFzLmhlaWdodCA9IHJlY3QuaGVpZ2h0O1xuXHRcdGNhbnZhcy5zdHlsZS5iYWNrZ3JvdW5kQ29sb3IgPSBcInJnYmEoMjUwLCAyNTAsIDI1MCwgMClcIjsgXG5cdFx0Y29udGFpbmVyLmFwcGVuZENoaWxkKGNhbnZhcyk7XG5cdH1cblx0ZWxzZXtcblx0XHR2YXIgY2FudmFzID0gY2FudmFzTGlzdFswXTtcblx0XHRjb25zb2xlLmxvZygnQ0FOVkFTJyk7XG5cdH1cblx0XG5cdHZhciBjb250ZXh0ID0gY2FudmFzLmdldENvbnRleHQoXCIyZFwiKTtcblx0Y29udGV4dC5jbGVhclJlY3QgKCAwICwgMCAsIGNhbnZhcy53aWR0aCwgY2FudmFzLmhlaWdodCApO1xuXHRcblxuXHRmdW5jdGlvbiB0aWNrKCkge1xuXHRcdHZhciB3ID0gY2FudmFzLndpZHRoO1xuXHRcdHZhciBoID0gY2FudmFzLmhlaWdodDtcblx0XHR2YXIgbW91c2UgPSBbbGFzdE1vdXNlTW92ZUV2ZW50LmNsaWVudFgvdywgbGFzdE1vdXNlTW92ZUV2ZW50LmNsaWVudFkvaF07XG5cdFx0Y29udGV4dC5zZXRUcmFuc2Zvcm0odywgMCwgMCwgaCwgMCwgMCk7XG5cdFx0Y29udGV4dC5maWxsU3R5bGUgPSBcInJnYmEoMjUwLCAyNTAsIDI1MCwgMC40KVwiO1xuXHRcdGNvbnRleHQuZmlsbFJlY3QoMCwwLHcsaCk7XG5cblx0XHQvLyBlZGdlc1xuXHRcdC8vIGNvbnRleHQuc3Ryb2tlU3R5bGUgPSBcIiMwMDBcIjtcblx0XHQvLyBmb3IodmFyIGk9MDsgaSA8IGVkZ2VzLmxlbmd0aDsgKytpKSB7XG5cdFx0Ly8gICAgIGNvbnRleHQubGluZVdpZHRoID0gMC4wMDAxO1xuXHRcdC8vICAgICB2YXIgZWRnZSA9IGVkZ2VzW2ldO1xuXHRcdC8vICAgICAvLyBpZiAoZWRnZS5waGVyb21vbiAhPSAwKXtcblx0XHQvLyAgICAgLy8gICAgIGNvbnRleHQubGluZVdpZHRoID0gTWF0aC5taW4oMC4wMDAwMSAqIGVkZ2UucGhlcm9tb24sIDAuMDEpO1xuXHRcdC8vICAgICAvLyB9IGVsc2Uge1xuXHRcdC8vICAgICAvLyAgICAgY29udGV4dC5saW5lV2lkdGggPSAwLjAwMDAxO1xuXHRcdC8vICAgICAvLyB9XG5cdFx0Ly8gICAgIGNvbnRleHQuYmVnaW5QYXRoKCk7XG5cdFx0Ly8gICAgIGNvbnRleHQubW92ZVRvKHBvaW50c1tlZGdlLnB0MS5pZF0ueCwgcG9pbnRzW2VkZ2UucHQxLmlkXS55KTtcblx0XHQvLyAgICAgY29udGV4dC5saW5lVG8ocG9pbnRzW2VkZ2UucHQyLmlkXS54LCBwb2ludHNbZWRnZS5wdDIuaWRdLnkpO1xuXHRcdC8vICAgICBjb250ZXh0LnN0cm9rZSgpO1xuXHRcdC8vIH1cblxuXHRcdC8vIC8vIHZlcnRpY2VzXG5cdFx0Ly8gZm9yKHZhciBpPTA7IGk8cG9pbnRzLmxlbmd0aDsgKytpKSB7XG5cdFx0Ly8gICAgIGNvbnRleHQuYmVnaW5QYXRoKClcblx0XHQvLyAgICAgdmFyIHBvaW50ID0gcG9pbnRzW2ldO1xuXHRcdC8vICAgICBpZiAoY2l0eVNldC5oYXMocG9pbnQuaWQpKSB7XG5cdFx0Ly8gICAgICAgICBjb250ZXh0LmZpbGxTdHlsZSA9IFwiIzAxMDFERlwiO1xuXHRcdC8vICAgICAgICAgY29udGV4dC5hcmMocG9pbnQueCwgcG9pbnQueSwgMC4wMDYsIDAsIDIqTWF0aC5QSSk7XG5cdFx0Ly8gICAgIH1cblx0XHQvLyAgICAgZWxzZSB7XG5cdFx0Ly8gICAgICAgICBjb250ZXh0LmZpbGxTdHlsZSA9IFwiIzAwMFwiO1xuXHRcdC8vICAgICAgICAgY29udGV4dC5hcmMocG9pbnRzW2ldLngsIHBvaW50c1tpXS55LCAwLjAwMywgMCwgMipNYXRoLlBJKTtcblx0XHQvLyAgICAgfVxuXHRcdC8vICAgICBjb250ZXh0LmNsb3NlUGF0aCgpO1xuXHRcdC8vICAgICBjb250ZXh0LmZpbGwoKTtcblx0XHQvLyB9XG5cblx0XHQvLyBtb3ZlIGFudHNcblx0XHRwb3B1bGF0aW9uLmZvckVhY2goZnVuY3Rpb24oYW50KXtcblx0XHRcdGFudC50cmFuc2l0KCk7XG5cdFx0fSk7XG5cblx0XHQvLyBmb3IgKGkgPSAwOyBpIDwgbmJBbnRzOyBpKyspIHtcblx0XHQvLyAgICAgcG9wdWxhdGlvbltpXS50cmFuc2l0KCk7XG5cdFx0Ly8gfVxuXG5cdFx0Ly8gcGhlcm9tb24gZXZhcG9yYXRpb25cblx0XHRlZGdlcy5mb3JFYWNoKGZ1bmN0aW9uKGVkZ2Upe1xuXHRcdFx0aWYoZWRnZS5waGVyb21vbiA+IDApe1xuXHRcdFx0XHRlZGdlLnBoZXJvbW9uIC09IDAuMDAwMTtcblx0XHRcdH1cblx0XHR9KTtcblxuXHRcdC8vIGZvciAoaSA9IDA7IGkgPCBlZGdlcy5sZW5ndGg7IGkrKykge1xuXHRcdFx0XG5cblx0XHQvLyBhbnRzXG5cdFx0cG9wdWxhdGlvbi5mb3JFYWNoKGZ1bmN0aW9uKGFudCl7XG5cdFx0XHRjb250ZXh0LmJlZ2luUGF0aCgpXG5cdFx0XHR2YXIgeCA9IGFudC54ICsgUkFORE9NTVZUKnJhbmRvbSgpO1xuXHRcdFx0dmFyIHkgPSBhbnQueSArIFJBTkRPTU1WVCpyYW5kb20oKTtcblxuXHRcdFx0Y29udGV4dC5maWxsU3R5bGUgPSBcImJsYWNrXCJcblx0XHRcdGNvbnRleHQuZmlsbFJlY3QoeCwgeSwgQU5UU0laRSwgQU5UU0laRSk7XG5cdFx0XHRjb250ZXh0LmNsb3NlUGF0aCgpO1xuXHRcdFx0Y29udGV4dC5maWxsKCk7XG5cdFx0fSlcblx0XHQvLyBmb3IodmFyIGk9MDsgaTxwb3B1bGF0aW9uLmxlbmd0aDsgKytpKSB7XG5cdFx0XHRcblx0XHQvLyB9XG5cblx0fTtcblx0XG5cdHZhciBsYXN0TW91c2VNb3ZlRXZlbnQgPSB7XG5cdFx0Y2xpZW50WDogMCxcblx0XHRjbGllbnRZOiAwXG5cdH07XG5cdFxuXHRjb250YWluZXIuYWRkRXZlbnRMaXN0ZW5lcignbW91c2Vtb3ZlJywgZnVuY3Rpb24oZSl7XG5cdFx0bGFzdE1vdXNlTW92ZUV2ZW50ID0gZTtcblx0fSk7XG5cdFxuXHR2YXIgcGF1c2VkID0gZmFsc2U7XG5cdFxuXHRmdW5jdGlvbiB0b2dnbGVQbGF5UGF1c2UoKXtcblx0XHRwYXVzZWQgPSAhcGF1c2VkO1xuXHRcdGlmKCFwYXVzZWQpXG5cdFx0XHRhbmltYXRlKCk7XG5cdH1cblxuXHRmdW5jdGlvbiByZXNldCgpe1xuXHRcdHBvcHVsYXRpb24gPSBbXTtcblx0XHRlZGdlcyA9IFtdO1xuXHRcdHBvaW50c0luZm9zID0gW107XG5cblx0XHRjYW5jZWxBbmltYXRpb25GcmFtZShhbmltSUQpO1xuXHR9XG5cdFxuXHQvLyBjb250YWluZXIuYWRkRXZlbnRMaXN0ZW5lcignY2xpY2snLCB0b2dnbGVQbGF5UGF1c2UpO1xuXG5cdGZ1bmN0aW9uIGFuaW1hdGUoKXtcblx0XHR0aWNrKCk7XG5cdFx0XG5cdFx0aWYoIXBhdXNlZClcblx0XHRcdGFuaW1JRCA9IHJlcXVlc3RBbmltYXRpb25GcmFtZShhbmltYXRlKTtcblx0fVxuXHRhbmltYXRlKCk7XG5cblxuXHRmdW5jdGlvbiBzZXRBbnRDb3VudChvcHRzKXtcblxuXHRcdHZhciBwcmV2aW91c0NvdW50ID0gcG9wdWxhdGlvbi5sZW5ndGg7XG5cblx0XHRjb25zb2xlLmxvZygnbmJBbnRzIDonLCBvcHRzLm5iQW50cyk7XG5cdFx0Y29uc29sZS5sb2coJ3ByZXZpb3VzQ291bnQgOicsIHByZXZpb3VzQ291bnQpO1xuXG5cblx0XHRpZiAob3B0cy5uYkFudHMgPiBwcmV2aW91c0NvdW50KXtcblx0XHRcdGNvbnNvbGUubG9nKCdwb2ludHNJbmZvcyA6JywgcG9pbnRzSW5mb3MpO1xuXHRcdFx0dmFyIEFudCA9IGFudEZ1bmN0aW9uKGNvbnRhaW5lciwgcG9pbnRzSW5mb3MsIG9wdHMpO1xuXG5cdFx0XHRmb3IgKHZhciBpID0gMDsgaSA8IG9wdHMubmJBbnRzIC0gcHJldmlvdXNDb3VudDsgaSsrKSB7XG5cdFx0XHRcdGNvbnNvbGUubG9nKCd0ZXN0Jyk7XG5cdFx0XHRcdHZhciBuZXdBbnQgPSBuZXcgQW50KEFudC5nZW5lcmF0ZVJhbmRTdGFydFBvaW50KCkpO1xuXHRcdFx0XHRuZXdBbnQuc2V0RGlyZWN0aW9uKCk7XG5cdFx0XHRcdHBvcHVsYXRpb24ucHVzaChuZXdBbnQpO1xuXHRcdFx0fVxuXHRcdH1cblx0XHRlbHNle1xuXHRcdFx0cG9wdWxhdGlvbiA9IHBvcHVsYXRpb24uc2xpY2UoMCwgb3B0cy5uYkFudHMpO1xuXHRcdFx0bmJBbnRzID0gcG9wdWxhdGlvbi5sZW5ndGg7XG5cdFx0XHRjb25zb2xlLmxvZygnTmIgQW50cyA6JywgcG9wdWxhdGlvbi5sZW5ndGgpO1xuXHRcdH1cblx0XHRjb25zb2xlLmxvZygnTm91dmVhdSBjb21wdGUgOicsIHBvcHVsYXRpb24ubGVuZ3RoKTtcblx0fVxuXG5cdGZ1bmN0aW9uIG1vZGlmeUFudHMob3B0cyl7XG5cdFx0c2V0QW50Q291bnQob3B0cyk7XG5cblx0XHRwb3B1bGF0aW9uLmZvckVhY2goZnVuY3Rpb24oYW50KXtcblx0XHRcdGFudC52ZWxvY2l0eSA9IG9wdHMudmVsb2NpdHk7XG5cdFx0XHRhbnQud2VpZ2h0ID0gb3B0cy53ZWlnaHQ7XG5cdFx0XHRhbnQucmVwU2l6ZSA9IG9wdHMucmVwU2l6ZTtcblx0XHRcdGFudC5yZXBTcGVlZCA9IG9wdHMucmVwU3BlZWQ7XG5cdFx0fSk7XG5cdH1cblx0XG5cdHJldHVybiB7XG5cdFx0dG9nZ2xlUGxheVBhdXNlOiB0b2dnbGVQbGF5UGF1c2UsXG5cdFx0cmVzZXQ6IHJlc2V0LFxuXHRcdC8vIHNob3VsZCBiZSBhIGdldHRlci9zZXR0ZXIsIGJ1dCBJRThcblx0XHRnZXRBbnRDb3VudDogZnVuY3Rpb24oKXtcblx0XHRcdHJldHVybiBwb3B1bGF0aW9uLmxlbmd0aDtcblx0XHR9LFxuXHRcdHNldEFudENvdW50OiBzZXRBbnRDb3VudCxcblx0XHRtb2RpZnlBbnRzOiBtb2RpZnlBbnRzXG5cdH1cbn1cbiIsIid1c2Ugc3RyaWN0JztcblxudmFyIHNxcnQgPSBNYXRoLnNxcnQ7XG52YXIgcG93ID0gTWF0aC5wb3c7XG5cbmZ1bmN0aW9uIHNpZ24oeCkge1xuXHRyZXR1cm4geCA/IHggPCAwID8gLTEgOiAxIDogMDtcbn1cblxuZnVuY3Rpb24gcmFuZ2Uoc3RhcnQsIGNvdW50KSB7XG4gICAgcmV0dXJuIEFycmF5LmFwcGx5KDAsIEFycmF5KGNvdW50KSkubWFwKGZ1bmN0aW9uIChlbGVtZW50LCBpbmRleCkge1xuICAgIFx0cmV0dXJuIGluZGV4ICsgc3RhcnRcbiAgICB9KTtcbn1cblxuZnVuY3Rpb24gZGlzdGFuY2UoYSwgYil7XG5cdHJldHVybiBzcXJ0KHBvdyhhLnggLSBiLngsIDIpICsgcG93KGEueSAtIGIueSwgMikpO1xufVxuXG5mdW5jdGlvbiBub3JtKHYpe1xuXHRyZXR1cm4gc3FydChwb3codi54LCAyKSArIHBvdyh2LnksIDIpKTtcbn1cblxubW9kdWxlLmV4cG9ydHMgPSB7XG5cdHNpZ246IHNpZ24sXG5cdHJhbmdlOiByYW5nZSxcblx0ZGlzdGFuY2U6IGRpc3RhbmNlLFxuXHRub3JtOiBub3JtXG59IiwiJ3VzZSBzdHJpY3QnXG5cbmZ1bmN0aW9uIFZlY3Rvcih4LCB5KSB7XG4gICAgdGhpcy54ID0geDsgICAgICAgICAgICAgICAgXG4gICAgdGhpcy55ID0geTtcbn1cblxuVmVjdG9yLnByb3RvdHlwZS5ub3JtID0gZnVuY3Rpb24oKXtcblx0cmV0dXJuIE1hdGguc3FydCh0aGlzLnggKiB0aGlzLnggKyB0aGlzLnkgKiB0aGlzLnkpO1xufVxuXG5WZWN0b3IucHJvdG90eXBlLm5vcm1hbGl6ZSA9IGZ1bmN0aW9uKCl7XG5cdHZhciBub3JtID0gdGhpcy5ub3JtKCk7XG5cdHRoaXMueCA9IHRoaXMueCAvIG5vcm07XG5cdHRoaXMueSA9IHRoaXMueSAvIG5vcm07XG59XG5cblxuXG5tb2R1bGUuZXhwb3J0cyA9IFZlY3RvcjsiXX0=
