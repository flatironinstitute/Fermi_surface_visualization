
/*

JQuery plugin populating 2 and 3 dimensional contours.

Requires nd_frame to be loaded.

Structure follows: https://learn.jquery.com/plugins/basic-plugin-creation/

*/
"use strict";

(function($) {

    var std_sizes_declarations = `
    int voxel_index;
    int i_block_num;
    int i_depth_num;
    int i_row_num;
    int i_col_num;
    float f_col_num;
    float f_row_num;
    float f_depth_num;
    vec3 location_offset;
    `;

    var get_sizes_macro = function(index_variable_name) {
        return `
        voxel_index = ${index_variable_name};
        
        // size of layer of rows and columns in 3d grid block
        int layer_voxels = uRowSize * uColSize;
        //int i_block_num;
        int block_index;

        if (uLayerSize > 1) {
            // possibly multiple grids in blocks.
            // size of block of rows/columns/layers
            int block_voxels = layer_voxels * uLayerSize;

            // block number for this voxel
            i_block_num = voxel_index / block_voxels;
            // ravelled index in block
            block_index = voxel_index - (i_block_num * block_voxels);
        } else {
            // only one block
            i_block_num = 0;
            block_index = voxel_index;
        }

        // instance depth of this layer
        i_depth_num = block_index / layer_voxels;
        // ravelled index in layer
        int i_layer_index = block_index - (i_depth_num * layer_voxels);

        i_row_num = i_layer_index/ uRowSize;
        i_col_num = i_layer_index - (i_row_num * uRowSize);

        f_col_num = float(i_col_num);
        f_row_num = float(i_row_num);
        f_depth_num = float(i_depth_num);
        //float f_block_num = float(i_block_num);  // not needed?
        location_offset = vec3(f_depth_num, f_row_num, f_col_num);
        `;
    };

    // functions to compute (x,y,z) location of offset relative to voxel location.
    var grid_location_decl = `
    vec3 grid_location(in vec3 offset) {
        vec3 rescaled = rescale_offset(offset);
        return grid_xyz(rescaled);
    }
    `;

    var locate_std_decl = `
    vec3 rescale_offset(in vec3 offset) {
        // convert voxel offset to block grid
        return location_offset + offset;
    }

    vec3 grid_xyz(in vec3 offset) {
        // convert block grid coords to xyz (trivial here)
        return offset;
    }
    ${grid_location_decl}
    `;

    // xxxx the scaling and and polar conversion could be separated eventually if useful.
    var locate_polar_scaled_decl = `

    // all samplers hold value in R component only.
    // [block, row] --> row_scaled
    uniform sampler2D RowScale;

    // [block, col] --> col_scaled
    uniform sampler2D ColumnScale;

    // [block, layer] --> layer_scaled
    uniform sampler2D LayerScale;

    float rescale_f(in float offset, in int index, in sampler2D scaling) {
        // note: indices are inverted from matrix notation matrix[y,x] === sampler(x,y) (???)
        //float x0 = texelFetch(scaling, ivec2(i_block_num, index), 0).r;
        //float x1 = texelFetch(scaling, ivec2(i_block_num, index+1), 0).r;
        float x0 = texelFetch(scaling, ivec2(index, i_block_num), 0).r;
        float x1 = texelFetch(scaling, ivec2(index+1, i_block_num), 0).r;
        return (x0 * (1.0 - offset)) + (x1 * offset);  // no clamping?
    }

    vec3 rescale_offset(in vec3 offset) {
        // convert voxel offset to block grid.
        // spherical coordinates using the "3rd major convention"
        // https://en.wikipedia.org/wiki/Spherical_coordinate_system#Conventions
        float r = rescale_f(offset[0], i_depth_num, LayerScale);
        // swapping phi and theta.
        float phi = rescale_f(offset[1], i_row_num, RowScale);
        float theta = rescale_f(offset[2], i_col_num, ColumnScale);
        return vec3(r, phi, theta);
    }

    vec3 grid_xyz(in vec3 spherical) {
        // convert block grid coords to xyz (trivial here)
        float r = spherical[0];
        // swapping phi and theta.
        float phi = spherical[1];
        float theta = spherical[2];
        //return vec3(r, theta, phi);
        
        float sint = sin(theta);
        float cost = cos(theta);
        float sinp = sin(phi);
        float cosp = cos(phi);
        float x = r * sinp * cost;
        float y = r * sinp * sint;
        float z = r * cosp;
        return vec3(x, y, z);
    }
    ${grid_location_decl}
    `;

    $.fn.webGL2crossingVoxels = function(options) {

        class WebGL2CrossingVoxels {
            constructor(options) {
                this.settings = $.extend({
                    feedbackContext: null,    // the underlying FeedbackContext context to use
                    valuesArray: null,   // the array buffer of values to contour
                    num_rows: null,
                    num_cols: null,
                    num_layers: 1,  // default to "flat"
                    num_blocks: 1,  // for physics simulations data may come in multiple blocks
                    grid_min: [0, 0, 0],
                    grid_max: [-1, -1, -1],  // disabled grid coordinate filtering (invalid limits)
                    rasterize: false,
                    threshold: 0,  // value at contour
                    // when getting compact arrays
                    // shrink the array sizes by this factor.
                    shrink_factor: 0.2,
                    location: "std",
                    // samplers are prepared by caller if needed.  Descriptors provided by caller.
                    samplers: {},
                    location_fill: -1e12,
                }, options);

                var s = this.settings;
                this.feedbackContext = s.feedbackContext;
                var nvalues = s.valuesArray.length;
                var nvoxels = s.num_rows * s.num_cols * s.num_layers * s.num_blocks;
                if (nvalues != nvoxels) {
                    // for now strict checking
                    throw new Error("voxels " + nvoxels + " don't match values " + nvalues);
                }
                // allocate and load buffer with a fresh name
                this.buffer = this.feedbackContext.buffer()
                this.buffer.initialize_from_array(s.valuesArray);
                var buffername = this.buffer.name;

                var vertex_shader;
                if (s.location == "std") {
                    vertex_shader = crossingVoxelsShader(locate_std_decl);
                } else if (s.location="polar_scaled") {
                    vertex_shader = crossingVoxelsShader(locate_polar_scaled_decl);
                } else {
                    throw new Error("unknown grid location type: " + s.location);
                }

                this.program = this.feedbackContext.program({
                    vertex_shader: vertex_shader,
                    fragment_shader: this.settings.fragment_shader,
                    feedbacks: {
                        index: {type: "int"},
                        location: {num_components: 3},
                        front_corners: {num_components: 4},
                        back_corners: {num_components: 4},
                    },
                });

                // set up input parameters
                //  indexing is [ix, iy, iz] -- z is fastest
                //var x_offset = 1;
                var z_offset = 1;
                var y_offset = s.num_cols;
                //var z_offset = s.num_cols * s.num_rows;
                var x_offset = s.num_cols * s.num_rows;
                var num_voxels = nvalues - (x_offset + y_offset + z_offset);

                var inputs = {};
                var add_input = function (ix, iy, iz) {
                    var name = (("a" + ix) + iy) + iz;
                    var dx = [0, x_offset][ix];
                    var dy = [0, y_offset][iy];
                    var dz = [0, z_offset][iz];
                    inputs[name] = {
                        per_vertex: true,
                        num_components: 1,
                        from_buffer: {
                            name: buffername,
                            skip_elements: dx + dy + dz,
                        }
                    }
                };
                add_input(0,0,0);
                add_input(0,0,1);
                add_input(0,1,0);
                add_input(0,1,1);
                add_input(1,0,0);
                add_input(1,0,1);
                add_input(1,1,0);
                add_input(1,1,1);

                this.runner = this.program.runner({
                    num_instances: 1,
                    vertices_per_instance: num_voxels,
                    rasterize: s.rasterize,
                    uniforms: {
                        // number of rows
                        uRowSize: {
                            vtype: "1iv",
                            default_value: [s.num_cols],
                        },
                        // numver of columns
                        uColSize: {
                            vtype: "1iv",
                            default_value: [s.num_rows],
                        },
                        // number of layers
                        uLayerSize: {
                            vtype: "1iv",
                            default_value: [s.num_layers],
                        },
                        // threshold value
                        uValue: {
                            vtype: "1fv",
                            default_value: [s.threshold],
                        },
                        u_grid_min: {
                            vtype: "3fv",
                            default_value: s.grid_min,
                        },
                        u_grid_max: {
                            vtype: "3fv",
                            default_value: s.grid_max,
                        },
                    },
                    inputs: inputs,
                    samplers: s.samplers,
                });
                this.front_corners_array = null;
                this.back_corners_array = null;
                this.index_array = null;
                this.compact_length = null;
            };
            run() {
                this.runner.install_uniforms();
                this.runner.run();
            };
            get_sphere_mesh(options) {
                // must be run after get_compacted_feedbacks has run at least once.
                var settings = $.extend({
                    THREE: null,   // required THREE instance
                    material: null, // material to use
                    radius: 1,  // shared radius for spheres
                    width_segments: 10,
                    height_segments: 10,}, options);
                settings.locations = this.compact_locations;
                return $.fn.webGL2crossingVoxels.spheresMesh(settings);
            };
            get_points_mesh(options) {
                // must be run after get_compacted_feedbacks has run at least once.
                var that = this;
                var settings = $.extend({
                    THREE: null,   // required THREE instance
                    size: null,
                    colorize: false,
                }, options);
                settings.locations = this.compact_locations;
                settings.center = this.compacted_feedbacks.mid;
                settings.radius = this.compacted_feedbacks.radius;
                if (settings.colorize) {
                    settings.colors = this.get_location_colors();
                }
                var result = $.fn.webGL2crossingVoxels.pointsMesh(settings);
                result.update_sphere_locations = function(locations, colors) {
                    locations = locations || that.compact_locations;
                    var geometry = result.geometry;
                    geometry.attributes.position.array = locations;
                    geometry.attributes.position.needsUpdate = true;
                    if (settings.colorize) {
                        colors = colors || that.get_location_colors();
                        geometry.attributes.color.array = colors;
                        geometry.attributes.color.needsUpdate = true;
                    }
                };
                return result;
            };
            get_compacted_feedbacks(location_only) {
                this.run();
                var location_fill = this.settings.location_fill;
                var rn = this.runner;
                if (!location_only) {
                    this.front_corners_array = rn.feedback_array(
                        "front_corners",
                        this.front_corners_array,
                    );
                    this.back_corners_array = rn.feedback_array(
                        "back_corners",
                        this.back_corners_array,
                    );
                    // xxxx locations are not always needed -- could optimize.
                    this.location_array = rn.feedback_array(
                        "location",
                        this.location_array,
                    );
                } else {
                    this.location_array = rn.feedback_array(
                        "location",
                        this.location_array,
                    );
                }
                this.index_array = rn.feedback_array(
                    "index",
                    this.index_array,
                );
                if (this.compact_length === null) {
                    // allocate arrays to size limit
                    this.compact_length = Math.floor(
                        this.settings.shrink_factor * this.index_array.length
                    );
                    this.compact_front_corners = new Float32Array(4 * this.compact_length);
                    this.compact_back_corners = new Float32Array(4 * this.compact_length);
                    this.compact_indices = new Int32Array(this.compact_length);
                    this.compact_locations = new Float32Array(3 * this.compact_length);
                }
                // compact the arrays
                this.compact_indices = this.feedbackContext.filter_degenerate_entries(
                    this.index_array, this.index_array, this.compact_indices, 1, -1
                );
                if (!location_only) {
                    this.compact_front_corners = this.feedbackContext.filter_degenerate_entries(
                        this.index_array, this.front_corners_array, this.compact_front_corners, 4, -1
                    );
                    this.compact_back_corners = this.feedbackContext.filter_degenerate_entries(
                        this.index_array, this.back_corners_array, this.compact_back_corners, 4, -1
                    );
                    // xxxx locations are not always needed -- could optimize.
                    this.compact_locations = this.feedbackContext.filter_degenerate_entries(
                        this.index_array, this.location_array, this.compact_locations, 3, location_fill
                    );
                } else {
                    this.compact_locations = this.feedbackContext.filter_degenerate_entries(
                        this.index_array, this.location_array, this.compact_locations, 3, location_fill
                    );
                }
                var mins = null;
                var maxes = null;
                var locs = this.compact_locations;
                var indices = this.compact_indices;
                if ((indices.length>0) && (indices[0]>=0)) {
                    mins = [locs[0], locs[1], locs[2]];
                    maxes = [locs[0], locs[1], locs[2]];
                    for (var i=0; i<indices.length; i++) {
                        if (indices[i]<0) {
                            break;
                        }
                        for (var k=0; k<3; k++) {
                            var v = locs[i*3 + k];
                            mins[k] = Math.min(mins[k], v);
                            maxes[k] = Math.max(maxes[k], v);
                        }
                    }
                }
                var n2 = 0;
                var mid = [];
                if (mins) {
                    for (var k=0; k<3; k++) {
                        mid.push(0.5 * (mins[k] + maxes[k]));
                        n2 += (mins[k] - maxes[k]) ** 2;
                    }
                }
                this.compacted_feedbacks = {
                    mid: mid,
                    radius: 0.5 * Math.sqrt(n2),
                    mins: mins,
                    maxes: maxes,
                    indices: this.compact_indices, 
                    front_corners: this.compact_front_corners,
                    back_corners: this.compact_back_corners,
                    locations: this.compact_locations,
                };
                return this.compacted_feedbacks;
            };
            get_location_colors() {
                var indices = this.compact_indices;
                var locations = this.compact_locations;
                var feedbacks = this.compacted_feedbacks;
                var mins = feedbacks.mins;
                var maxes = feedbacks.maxes;
                var colors = this.compact_colors;
                if (!colors) {
                    colors = new Float32Array(locations.length);
                    this.compact_colors = colors;
                }
                if ((!indices) || (indices[0] < 0)) {
                    return colors;  // no points: do nothing
                }
                var diffs = [];
                var base_intensity = 0.2;
                for (var j=0; j<3; j++) {
                    var d = maxes[j] - mins[j];
                    if (d < 1e-9) {
                        d = 1.0;
                    }
                    diffs.push(d / (1 - base_intensity));
                }
                for (var i=0; i<indices.length; i++) {
                    if (indices[i]<0) {
                        break;
                    }
                    for (var j=0; j<3; j++) {
                        var ij = i * 3 + j;
                        colors[ij] = base_intensity + (locations[ij] - mins[j])/diffs[j];
                    }
                }
                return colors;
            };
            reset_three_camera(camera, radius_multiple, orbit_control) {
                // adjust three.js camera to look at current voxels
                var cf = this.compacted_feedbacks;
                if ((!cf) || (!cf.mins)) {
                    // no points -- punt
                    return;
                }
                var cx = cf.mid[0];
                var cy = cf.mid[1];
                var cz = cf.mid[2];
                radius_multiple = radius_multiple || 3;
                camera.position.x = cx;
                camera.position.y = cy;
                camera.position.z = cz + radius_multiple * cf.radius;
                camera.lookAt(cf.mid[0], cf.mid[1], cf.mid[2]);
                if (orbit_control) {
                    orbit_control.center.x = cx;
                    orbit_control.center.y = cy;
                    orbit_control.center.z = cz;
                }
                return camera;
            }
            set_threshold(value) {
                this.runner.change_uniform("uValue", [value]);
            };
            set_grid_limits(grid_mins, grid_maxes) {
                this.runner.change_uniform("u_grid_min", grid_mins);
                this.runner.change_uniform("u_grid_max", grid_maxes);
            };
        };

        var crossingVoxelsShader = function(grid_location_declaration) {
            return `#version 300 es

        // global length of rows
        uniform int uRowSize;

        // global number of columnss
        uniform int uColSize;

        // global number of layers (if values are in multiple blocks, else 0)
        uniform int uLayerSize;
        
        // global contour threshold
        uniform float uValue;

        // global grid thresholds
        //  (I tried integers but it didn't work, couldn't debug...)
        uniform vec3 u_grid_min, u_grid_max;

        // per mesh function values at voxel corners
        in float a000, a001, a010, a011, a100, a101, a110, a111;

        // corners feedbacks
        out vec4 front_corners, back_corners;

        // location feedback
        out vec3 location;

        // index feedback
        flat out int index;

        ${std_sizes_declarations}
        ${grid_location_declaration}

        void main() {
            // default to invalid index indicating the voxel does not cross the value.
            index = -1;
            front_corners = vec4(a000, a001, a010, a011);
            back_corners = vec4(a100, a101, a110, a111);

            ${get_sizes_macro("gl_VertexID")}
            //location = location_offset;
            vec3 rescaled = rescale_offset(vec3(0,0,0));
            //location = grid_location(vec3(0,0,0));
            location = grid_xyz(rescaled);

            bool voxel_ok = true;
            if (u_grid_min[0] < u_grid_max[0]) {
                // voxel coordinate filtering is enabled
                voxel_ok = ( 
                    (u_grid_min[0] <= rescaled[0]) && (rescaled[0] < u_grid_max[0]) &&
                    (u_grid_min[1] <= rescaled[1]) && (rescaled[1] < u_grid_max[1]) &&
                    (u_grid_min[2] <= rescaled[2]) && (rescaled[2] < u_grid_max[2]) );
            }

            // Dont tile last column/row/layer which wraps around
            if ((voxel_ok) && 
                (i_col_num < (uRowSize - 1)) && 
                (i_row_num < (uColSize - 1)) &&
                (i_depth_num < (uLayerSize - 1))) {
                float m = front_corners[0];
                float M = front_corners[0];
                vec4 corners = front_corners;
                for (int j=0; j<2; j++) {
                    for (int i=0; i<4; i++) {
                        float c = corners[i];
                        m = min(c, m);
                        M = max(c, M);
                    }
                    corners = back_corners;
                }
                if ((m <= uValue) && (M > uValue)) {
                    // the pixel crosses the threshold
                    index = gl_VertexID;
                }
            }
        }
        `;};
        return new WebGL2CrossingVoxels(options);
    };

    $.fn.webGL2crossingVoxels.pointsMesh = function (options) {
        var settings = $.extend({
            THREE: null,   // required THREE instance
            locations: null,  // inifial points locations, required
            colors: null, // optional
            radius: 1.0,  // radius of bounding sphere
            center: [0, 0, 0],  // center of bounding sphere
            size: null,
        }, options);
        var THREE = settings.THREE;
        var locations = settings.locations;
        var c = settings.center;
        var size = settings.size || settings.radius * 0.01;
        var geometry = new THREE.BufferGeometry();
        geometry.setAttribute( 'position', new THREE.Float32BufferAttribute( locations, 3 ) );
        geometry.boundingSphere = new THREE.Sphere(new THREE.Vector3(c[0], c[1], c[2]), settings.radius);
        var vertex_colors = false;
        if (settings.colors) {
            geometry.setAttribute( 'color', new THREE.Float32BufferAttribute( settings.colors, 3 ) );
            vertex_colors = true;
        }
        var material = new THREE.PointsMaterial( { size: size, vertexColors: vertex_colors } );
        var points = new THREE.Points( geometry, material );
        points.update_sphere_locations = function(locations) {
            geometry.attributes.position.array = locations;
            geometry.attributes.position.needsUpdate = true;
        };
        return points;
    };

    $.fn.webGL2crossingVoxels.spheresMesh = function (options) {
        var settings = $.extend({
            THREE: null,   // required THREE instance
            material: null, // material to use
            locations: null,  // inifial sphere locations
            radius: 1,  // shared radius for spheres
            width_segments: 10,
            height_segments: 10,
        }, options);
        var THREE = settings.THREE;
        var locations = settings.locations;
        var geometry = new THREE.SphereBufferGeometry( settings.radius, settings.width_segments, settings.height_segments);
        var count = Math.floor(locations.length/3);
        var mesh = new THREE.InstancedMesh( geometry, settings.material, count );
        mesh.update_sphere_locations = function(locations) {
            var matrixArray = mesh.instanceMatrix.array;
            var translation_offset = 12;
            var matrix_size = 16;
            for (var i=0; i<count; i++) {
                var matrixStart = i * matrix_size + translation_offset;
                var locationStart = i * 3;
                // copy the translation portion of the matrix from the location.
                for (var j=0; j<3; j++) {
                    matrixArray[matrixStart + j] = locations[locationStart + j];
                }
            }
            mesh.instanceMatrix.needsUpdate = true;
        };
        // set up all matrices
        var M = new THREE.Matrix4();
        for (var i=0; i<count; i++) {
            mesh.setMatrixAt( i, M );
        }
        mesh.update_sphere_locations(locations);
        return mesh;
    };

    $.fn.webGL2crossingVoxels.example = function (container) {
        var gl = $.fn.feedWebGL2.setup_gl_for_example(container);

        var context = container.feedWebGL2({
            gl: gl,
        });
        var valuesArray = new Float32Array([
            0,0,0,
            0,0,0,
            0,0,0,

            0,0,0,
            0,1,0,
            0,0,0,

            0,0,0,
            0,0,0,
            0,0,0,
        ]);
        var crossing = container.webGL2crossingVoxels({
            feedbackContext: context,
            valuesArray: valuesArray,
            num_rows: 3,
            num_cols: 3,
            num_layers: 3, 
            threshold: 0.5,
            shrink_factor: 0.8,
        });
        var compacted = crossing.get_compacted_feedbacks();
        var indices = compacted.indices;
        var front_corners = compacted.front_corners;
        var back_corners = compacted.back_corners;
        var ci = 0;
        for (var i=0; i<indices.length; i++) {
            $("<br/>").appendTo(container);
            $("<span> " + indices[i] + " </span>").appendTo(container);
            for (var j=0; j<4; j++) {
                $("<span> " + front_corners[ci] + " " + back_corners[ci] + " </span>").appendTo(container);
                ci ++;
            }
        }
    };

    $.fn.webGL2TriangulateVoxels = function (options) {
        class WebGL2TriangulateVoxels {

            constructor(options) { 
                this.settings = $.extend({
                    feedbackContext: null,
                    // array of indices (from crossing voxels)
                    indices: null,
                    // array of corners (from crossing voxels)
                    front_corners: null,
                    back_corners: null,
                    // volume dimensions
                    num_rows: null,
                    num_cols: null,
                    num_layers: 0,  // if >1 then indexing in multiple blocks
                    dx: [1, 0, 0],
                    dy: [0, 1, 0],
                    dz: [0, 0, 1],
                    translation: [0, 0, 0],
                    color: [1, 1, 1],
                    rasterize: false,
                    threshold: 0,  // value at contour
                    invalid_coordinate: -100000,  // invalidity marker for positions
                    location: "std",
                    // samplers are prepared by caller if needed.  Descriptors provided by caller.
                    samplers: {},
                }, options);
                var s = this.settings;
                this.feedbackContext = s.feedbackContext;

                // allocate and load buffers with a fresh name
                this.index_buffer = this.feedbackContext.buffer()
                this.index_buffer.initialize_from_array(s.indices);
                this.front_corner_buffer = this.feedbackContext.buffer()
                this.front_corner_buffer.initialize_from_array(s.front_corners);
                this.back_corner_buffer = this.feedbackContext.buffer()
                this.back_corner_buffer.initialize_from_array(s.back_corners);
                // add vertex count bogus input for Firefox
                const N_TETRAHEDRA = 6;
                const N_TRIANGLES = 2;  
                const N_VERTICES = 3;
                var vertices_per_instance = N_TETRAHEDRA * N_TRIANGLES * N_VERTICES;
                this.vertices_per_instance = vertices_per_instance;
                // add vertex count bogus input for Firefox
                var vertexNumArray = new Float32Array(Array.from(Array(vertices_per_instance).keys()));
                this.vertex_num_buffer = this.feedbackContext.buffer()
                this.vertex_num_buffer.initialize_from_array(vertexNumArray);

                var vertex_shader;
                if (s.location == "std") {
                    vertex_shader = triangulate_vertex_shader(locate_std_decl);
                } else if (s.location="polar_scaled") {
                    vertex_shader = triangulate_vertex_shader(locate_polar_scaled_decl);
                    //vertex_shader = triangulate_vertex_shader(locate_std_decl);
                } else {
                    throw new Error("unknown grid location type: " + s.location);
                }

                this.program = this.feedbackContext.program({
                    vertex_shader: vertex_shader,
                    fragment_shader: tetrahedra_fragment_shader,
                    feedbacks: {
                        vPosition: {num_components: 3},
                        vNormal: {num_components: 3},
                        vColor: {num_components: 3},
                    },
                })

                this.runner = this.program.runner({
                    run_type: "TRIANGLES",
                    num_instances: s.indices.length,
                    vertices_per_instance: vertices_per_instance,
                    rasterize: s.rasterize,
                    uniforms: {
                        uRowSize: {
                            vtype: "1iv",
                            default_value: [s.num_cols],
                        },
                        uColSize: {
                            vtype: "1iv",
                            default_value: [s.num_rows],
                        },
                        // number of layers
                        uLayerSize: {
                            vtype: "1iv",
                            default_value: [s.num_layers],
                        },
                        uValue: {
                            vtype: "1fv",
                            default_value: [s.threshold],
                        },
                        dx: {
                            vtype: "3fv",
                            default_value: s.dx,
                        },
                        dy: {
                            vtype: "3fv",
                            default_value: s.dy,
                        },
                        dz: {
                            vtype: "3fv",
                            default_value: s.dz,
                        },
                        translation: {
                            vtype: "3fv",
                            default_value: s.translation,
                        },
                        uInvalid: {
                            vtype: "1fv",
                            default_value: [s.invalid_coordinate],
                        },
                    },
                    inputs: {
                        index: {
                            per_vertex: false,
                            num_components: 1,
                            type: "int",
                            from_buffer: {
                                name: this.index_buffer.name,
                            },
                        },
                        front_corners: {
                            per_vertex: false,
                            num_components: 4,
                            from_buffer: {
                                name: this.front_corner_buffer.name,
                            },
                        },
                        back_corners: {
                            per_vertex: false,
                            num_components: 4,
                            from_buffer: {
                                name: this.back_corner_buffer.name,
                            },
                        },
                        aVertexCount: {   // bogus attribute required by Firefox
                            per_vertex: true,
                            num_components: 1,
                            from_buffer: {
                                name: this.vertex_num_buffer.name,
                            },
                        },
                    },
                    samplers: s.samplers,
                });
            };
            run() {
                this.runner.install_uniforms();
                this.runner.run();
            };
            set_threshold(value) {
                //this.runner.uniforms.uValue.value = [value];
                this.runner.change_uniform("uValue", [value]);
                //this.runner.run();
            };
            get_positions(optionalPreAllocatedArrBuffer) {
                return this.runner.feedback_array(
                    "vPosition",
                    optionalPreAllocatedArrBuffer);
            };
            get_normals(optionalPreAllocatedArrBuffer) {
                return this.runner.feedback_array(
                    "vNormal",
                    optionalPreAllocatedArrBuffer);
            };
            get_colors(optionalPreAllocatedArrBuffer) {
                return this.runner.feedback_array(
                    "vColor",
                    optionalPreAllocatedArrBuffer);
            };
        };

        var triangulate_vertex_shader = function(grid_location_declaration) {
            return `#version 300 es

        // global length of rows, cols inputs
        uniform int uRowSize;
        uniform int uColSize;
        // global number of layers (if values are in multiple blocks, else 0)
        uniform int uLayerSize;
        
        // global contour threshold input
        uniform float uValue;
        
        // uniform offsets in xyz directions
        uniform vec3 dx, dy, dz, translation;
        
        // invalid value marker
        uniform float uInvalid;

        // per mesh corner values
        in vec4 front_corners, back_corners;

        // per mesh ravelled voxel index
        in int index;

        // bogus vertex attribute required by Firefox (but not Chrome)
        in float aVertexCount;

        // feedbacks out
        out vec3 vColor, vPosition, vNormal;

        // Which vertex in which triangle on which tetrahedron?
        //   gl_VertexID encodes tetrahedron_number 0..4, triangle_number 0..1, vertex number 0..2
        //   for a total of 6 * 2 * 3 = 36 vertices per "mesh instance".
        //   aVertexCount = tetrahedron_number * 10 + triangle_number * 3 + vertex_number;
        const int N_TETRAHEDRA = 6; // tetrahedra per cube
        const int N_TRIANGLES = 2;  // triangles per tetrahedron
        const int N_VERTICES = 3;   // vertices per triangle
        const int N_CORNERS = 8;    // number of cube corners
        const int N_T_VERTICES = 4; // number of vertices in a tetrahedron

        // Crossing index is binary integer associated with each tetrahedron of form
        //   (triangle_num << 4) || ((fA > v) << 3 || ((fB > v) << 2 || ((fC > v) << 1 || ((fD > v)
        const int N_CROSSING_INDEX = 32;  // 2 ** 5

        // corner offsets
        const vec3 offsets[N_CORNERS] = vec3[] (
            vec3(0.0, 0.0, 0.0),
            vec3(0.0, 0.0, 1.0),
            vec3(0.0, 1.0, 0.0),
            vec3(0.0, 1.0, 1.0),
            vec3(1.0, 0.0, 0.0),
            vec3(1.0, 0.0, 1.0),
            vec3(1.0, 1.0, 0.0),
            vec3(1.0, 1.0, 1.0)
        );

        // vertex indices for tiling tetrahedra per tetrahedron number
        const int A_INDEX = 0;
        const int[N_TETRAHEDRA] B_index = int[] (4, 6, 2, 3, 1, 5) ;
        const int[N_TETRAHEDRA] C_index = int[] (6, 2, 3, 1, 5, 4) ;
        const int D_INDEX = 7;

        // crossing index to triangle vertices endpoint indices
        const int U_left[N_CROSSING_INDEX] = int[] (
            -1, 3, 2, 0, 1, 0, 0, 0, 0, 1, 1, 1, 2, 2, 3,-1,
            -1,-1,-1, 0,-1, 0, 0,-1,-1, 1, 1,-1, 2,-1,-1,-1);
        const int U_right[N_CROSSING_INDEX] = int[] (
            -1, 1, 1, 2, 2, 1, 1, 2, 2, 0, 0, 2, 0, 1, 1,-1,
            -1,-1,-1, 3,-1, 3, 2,-1,-1, 3, 2,-1, 1,-1,-1,-1);
        const int V_left[N_CROSSING_INDEX] = int[] (
            -1, 3, 2, 1, 1, 0, 3, 0, 0, 2, 1, 1, 3, 2, 3,-1,
            -1,-1,-1, 1,-1, 2, 3,-1,-1, 2, 3,-1, 3,-1,-1,-1);
        const int V_right[N_CROSSING_INDEX] = int[] (
            -1, 0, 3, 2, 0, 3, 1, 1, 3, 0, 2, 3, 0, 0, 2,-1,
            -1,-1,-1, 2,-1, 3, 1,-1,-1, 0, 2,-1, 0,-1,-1,-1);
        const int W_left[N_CROSSING_INDEX] = int[] (
            -1, 3, 2, 0, 1, 2, 0, 0, 0, 1, 3, 1, 2, 2, 3,-1,
            -1,-1,-1, 1,-1, 2, 3,-1,-1, 2, 3,-1, 3,-1,-1,-1);
        const int W_right[N_CROSSING_INDEX] = int[] (
            -1, 2, 0, 3, 3, 1, 2, 3, 1, 3, 0, 0, 1, 3, 0,-1,
            -1,-1,-1, 3,-1, 1, 2,-1,-1, 3, 0,-1, 1,-1,-1,-1);

        ${std_sizes_declarations}
        ${grid_location_declaration}

        void main() {

            // initially set output point to invalid
            gl_Position = vec4(uInvalid, uInvalid, uInvalid, uInvalid);
            vPosition = gl_Position.xyz;
            // use the bogus vertexCount parameter so it is not erased by the optimizer
            float grey = aVertexCount / float(N_TETRAHEDRA * N_TRIANGLES * N_VERTICES);
            vColor = vec3(float(gl_VertexID) * 0.01, grey, 0.0);  // temp value for debugging
            vNormal = vec3(0.0, 0.0, 1.0);    // arbitrary initial value

            ${get_sizes_macro("index")}

            // Dont tile last column which wraps around or last row
            if ((index >= 0) && (i_col_num < (uRowSize - 1)) && (i_row_num < (uColSize - 1))) {
                // determine which vertex in which triangle in which tetrahedron to interpolate
                int iVertexCount = gl_VertexID;
                int iTetrahedronNumber = iVertexCount / (N_TRIANGLES * N_VERTICES);
                int iTetIndex = iVertexCount - (N_TRIANGLES * N_VERTICES) * iTetrahedronNumber;
                int iTriangleNumber = iTetIndex / N_VERTICES;
                int iVertexNumber = iTetIndex - (iTriangleNumber * N_VERTICES);
                // offsets of vertices for this tet number
                vec3 t_offsets[N_T_VERTICES] = vec3[](
                    offsets[A_INDEX],
                    offsets[B_index[iTetrahedronNumber]],
                    offsets[C_index[iTetrahedronNumber]],
                    offsets[D_INDEX]
                );
                vec3[4] grid_locations = vec3[](
                    grid_location(t_offsets[0]),
                    grid_location(t_offsets[1]),
                    grid_location(t_offsets[2]),
                    grid_location(t_offsets[3])
                );
                // weights as array
                float wts[N_CORNERS] = float[](
                    front_corners[0], front_corners[1], front_corners[2], front_corners[3], 
                    back_corners[0], back_corners[1], back_corners[2], back_corners[3]);
                    // a000, a001, a010, a011, a100, a101, a110, a111
                // weights of vertices for this tet number
                float t_wts[N_T_VERTICES] = float[](
                    wts[A_INDEX],
                    wts[B_index[iTetrahedronNumber]],
                    wts[C_index[iTetrahedronNumber]],
                    wts[D_INDEX]
                );
                // crossing index
                int ci = iTriangleNumber << 1;
                if (t_wts[0] > uValue) { ci = ci + 1; }
                ci = ci << 1;
                if (t_wts[1] > uValue) { ci = ci + 1; }
                ci = ci << 1;
                if (t_wts[2] > uValue) { ci = ci + 1; }
                ci = ci << 1;
                if (t_wts[3] > uValue) { ci = ci + 1; }

                // If U_left[ci] for this corner is negative (invalid index) then there is no such triangle here.
                if (U_left[ci] >= 0) {
                    int SegLs[N_VERTICES] = int[](U_left[ci], V_left[ci], W_left[ci]);
                    int SegRs[N_VERTICES] = int[](U_right[ci], V_right[ci], W_right[ci]);
                    vec3[N_VERTICES] combined_offsets;
                    // compute intercepts for all vertices of the triangle
                    for (int vnum=0; vnum<N_VERTICES; vnum++) {
                        int SegL = SegLs[vnum];
                        int SegR = SegRs[vnum];
                        vec3 offsetL = grid_locations[SegL];
                        vec3 offsetR = grid_locations[SegR];
                        float wtL = t_wts[SegL];
                        float wtR = t_wts[SegR];
                        // check denominator is not too small? xxxx
                        float delta = (wtL - uValue) / (wtL - wtR);
                        combined_offsets[vnum] = ((1.0 - delta) * offsetL) + (delta * offsetR);
                    }
                    vec3 vertex = combined_offsets[iVertexNumber];
                    vPosition = dx * vertex[0] + dy * vertex[1] + dz * vertex[2] + translation;
                    gl_Position.xyz = vPosition;
                    gl_Position[3] = 1.0;
                    //vdump = float[4](vertex[0], vertex[1], vertex[2], delta);

                    vec3 nm = cross(combined_offsets[1] - combined_offsets[0], combined_offsets[2] - combined_offsets[0]);
                    float ln = length(nm);
                    if (ln > 1e-12) {
                        vNormal = nm / ln;
                    }
                    vColor = abs(vNormal);  // XXX FOR TESTING ONLY
                }
            }
            //vPosition = gl_Position.xyz;
        }
        `;};

        var tetrahedra_fragment_shader = `#version 300 es
        #ifdef GL_ES
            precision highp float;
        #endif
        in vec3 vColor;
        out vec4 color;

        void main() {
            color = vec4(vColor, 1.0);
        }
        `;

        return new WebGL2TriangulateVoxels(options);
    };

    $.fn.webGL2TriangulateVoxels.example = function (container) {
        var gl = $.fn.feedWebGL2.setup_gl_for_example(container);

        var context = container.feedWebGL2({
            gl: gl,
        });
        /*
        Copied from data dump:

        0 0 0 0 0 0 0 0 1
        1 0 0 0 0 0 0 1 0
        3 0 0 0 1 0 0 0 0
        4 0 0 1 0 0 0 0 0
        9 0 0 0 0 0 1 0 0
        10 0 0 0 0 1 0 0 0
        12 0 1 0 0 0 0 0 0
        13 1 0 0 0 0 0 0 0
        -1 -1 -1 -1 -1 -1 -1 -1 -1
        -1 -1 -1 -1 -1 -1 -1 -1 -1
        -1 -1 -1 -1 -1 -1 -1 -1 -1

        First column is index, alternating columns are front/back corners after.
        */
        var front_cornersArray = new Float32Array([
            0, 0, 0, 0,
            0, 0, 0 ,1,
            0, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 0,
            1, 0, 0, 0,
            -1.,-1,-1,-1, 
            -1.,-1,-1,-1, 
            -1.,-1,-1,-1, 
        ]);
        var back_cornersArray = new Float32Array([
            0, 0, 0, 1,
            0, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 0,
            1, 0, 0, 0,
            0, 0, 0, 0,
            -1.,-1,-1,-1, 
            -1.,-1,-1,-1, 
            -1.,-1,-1,-1, 
        ]);
        var indexArray = new Int32Array([
            0,1,3,4,9,10,12,13,-1,-1,-1
        ]);
        var segments = container.webGL2TriangulateVoxels({
            feedbackContext: context,
            indices: indexArray,
            front_corners: front_cornersArray,
            back_corners: back_cornersArray,
            num_rows: 3,
            num_cols: 3,
            rasterize: true,
            dx: [0.3, 0, 0],
            dy: [0, 0.3, 0],
            dz: [0, 0, 0.3],
            threshold: 0.5,
        });
        segments.run();
        var positions = segments.get_positions();
        for (var i=0; i<positions.length; i++) {
            if (i % 3 == 0) {
                $("<br/>").appendTo(container);
            }
            $("<span> " + positions[i] + " </span>").appendTo(container);
        }
    };

    $.fn.webGL2surfaces3dopt = function (options) {
        // "optimized surfaces" by truncating buffer sizes
        // which may result in some data omission in dense cases.
        class WebGL2Surfaces3dOpt {
            constructor(options) {
                var that = this;
                this.settings = $.extend({
                    // default settings:
                    shrink_factor: 0.1, // how much to shrink buffers
                    feedbackContext: null,    // the underlying FeedbackContext context to use
                    valuesArray: null,   // the array buffer of values to contour
                    num_rows: null,
                    num_cols: null,
                    num_layers: 1,  // default to "flat"
                    num_blocks: 1,
                    //dx: [1, 0, 0],
                    //dy: [0, 1, 0],
                    //dz: [0, 0, 1],
                    //translation: [-1, -1, 0],
                    color: [1, 1, 1],
                    rasterize: false,
                    threshold: 0,  // value at contour
                    invalid_coordinate: -100000,  // invalidity marker for positions, must be very negative
                    grid_min: [0, 0, 0],
                    grid_max: [-1, -1, -1],  // disabled grid coordinate filtering (invalid limits)
                    after_run_callback: null,   // call this after each run.
                    // method of conversion from grid coordinates to world coordinates
                    location: "std", 
                    // parameters needed by location method if any.
                    location_parameters: null,
                }, options);
                this.check_geometry();
                var s = this.settings;
                this.feedbackContext = s.feedbackContext;
                var container = $(this.feedbackContext.canvas);
                if (!this.feedbackContext) {
                    throw new Error("Feedback context required.");
                }
                var nvalues = s.valuesArray.length;
                var nvoxels = s.num_rows * s.num_cols * s.num_layers * s.num_blocks;
                if (nvalues != nvoxels) {
                    // for now strict checking
                    throw new Error("voxels " + nvoxels + " don't match values " + nvalues);
                }
                // samplers for location conversion, if any
                this.samplers = {};
                this.textures = {}
                if (s.location == "polar_scaled") {
                    // set up scaling textures
                    this.samplers.RowScale = this.feedbackContext.texture("RowScale", "FLOAT", "RED", "R32F");
                    var set_up_sampler = function(name, size) {
                        var texture = that.feedbackContext.texture(name, "FLOAT", "RED", "R32F");
                        texture.load_array(s.location_parameters[name], size, s.num_blocks)
                        that.textures[name] = texture;
                        that.samplers[name] = {dim: "2D", from_texture: name};
                    };
                    set_up_sampler("RowScale", s.num_rows+1);
                    set_up_sampler("ColumnScale", s.num_cols+1);
                    set_up_sampler("LayerScale", s.num_layers+1);
                }
                this.crossing = container.webGL2crossingVoxels({
                    feedbackContext: this.feedbackContext,
                    valuesArray: s.valuesArray,
                    num_rows: s.num_rows,
                    num_cols: s.num_cols,
                    num_layers: s.num_layers,
                    num_blocks: s.num_blocks,
                    threshold: s.threshold,
                    shrink_factor: s.shrink_factor,
                    grid_min: s.grid_min,
                    grid_max: s.grid_max,  // disabled grid coordinate filtering (invalid limits)
                    location: s.location,
                    samplers: this.samplers,
                    // never rasterize the crossing pixels
                });
                // initialize segmenter upon first run.
                this.segments = null; 
            };
            check_geometry() {
                // arrange the geometry parameters to fit in [-1:1] cube unless specified otherwise
                var s = this.settings;
                if (s.location != "std") {
                    return;  // don't mess with non-standard geometry
                }
                if (!s.dx) {
                    // geometry needs specifying:
                    var max_dimension = Math.max(s.num_rows, s.num_cols, s.num_layers);
                    var dpixel = 2.0 / max_dimension;
                    s.dx = [dpixel, 0, 0];
                    s.dy = [0, dpixel, 0];
                    s.dz = [0, 0, dpixel];
                    if (!s.translation) {
                        s.translation = [-0.5 * s.num_cols * dpixel, -0.5 * s.num_rows * dpixel, -0.5 * s.num_layers * dpixel]
                    }
                }
            }
            run () {
                var s = this.settings;
                var compacted = this.crossing.get_compacted_feedbacks();
                if (!this.segments) {
                    var container = $(this.feedbackContext.canvas);
                    this.segments = container.webGL2TriangulateVoxels({
                        feedbackContext: this.feedbackContext,
                        indices: compacted.indices,
                        front_corners: compacted.front_corners,
                        back_corners: compacted.back_corners,
                        num_rows: s.num_rows,
                        num_cols: s.num_cols,
                        num_layers: s.num_layers,
                        num_blocks: s.num_blocks,
                        rasterize: s.rasterize,
                        dx: s.dx,
                        dy: s.dy,
                        dz: s.dz,
                        translation: s.translation,
                        threshold: s.threshold,
                        invalid_coordinate: s.invalid_coordinate,
                        location: s.location,
                        samplers: this.samplers,
                    });
                } else {
                    // reset buffer content
                    this.segments.index_buffer.copy_from_array(
                        compacted.indices
                    );
                    this.segments.front_corner_buffer.copy_from_array(
                        compacted.front_corners
                    );
                    this.segments.back_corner_buffer.copy_from_array(
                        compacted.back_corners
                    );
                }
                this.indices = compacted.indices;
                this.vertices_per_instance = this.segments.vertices_per_instance;
                this.segments.run();
                var after_run_callback = this.settings.after_run_callback;
                if (after_run_callback) {
                    after_run_callback(this);
                }
                //var positions = segments.get_positions();
            };
            colorization(voxel_color_source, vertex_color_destination) {
                // apply voxel colors to vertices for active voxel indices
                var indices = this.indices;
                var vertices_per_instance = this.vertices_per_instance;
                var skip_index = vertices_per_instance * 3;
                var num_indices = indices.length;
                var count = 0;
                for (var i=0; i<num_indices; i++){
                    var index = indices[i];
                    if (index < 0) {
                        count += skip_index;
                    } else {
                        var cindex = 3 * index;
                        for (var vn=0; vn<vertices_per_instance; vn++) {
                            for (var cn=0; cn<3; cn++) {
                                vertex_color_destination[count] = voxel_color_source[cindex + cn];
                                count ++;
                            }
                        }
                    }
                }
                return vertex_color_destination;
            };
            linked_three_geometry (THREE, clean, normal_binning) {
                // create a three.js geometry linked to the current positions feedback array.
                // xxxx only one geometry may be linked at a time.
                // this is a bit convoluted in an attempt to only update attributes when needed.
                var that = this;
                var positions, normals;
                if (clean) {
                    var pn = this.clean_positions_and_normals(normal_binning);
                    positions = pn.positions;
                    normals = pn.normals;
                } else {
                    positions = this.get_positions();
                    normals = this.get_normals();
                }
                var colors = this.get_colors();  // xxxx remove this? (debug only)
                var geometry = new THREE.BufferGeometry();
                geometry.setAttribute( 'position', new THREE.BufferAttribute( positions, 3 ) );
                geometry.setAttribute( 'normal', new THREE.BufferAttribute( normals, 3 ) );
                geometry.setAttribute( 'color', new THREE.BufferAttribute( colors, 3 ) );
                that.link_needs_update = false;
                var after_run = function(that) {
                    debugger;
                    that.link_needs_update = true;
                }
                var check_update_link = function(nbins) {
                    // update the surface if needed, using nbins for normal_binning if provided.
                    var do_clean = clean || nbins;
                    var bin_size = nbins || normal_binning;
                    // update the geometry positions array in place and mark for update in geometry
                    if ((!that.link_needs_update) && (!nbins)) {
                        // only update upon request and only if needed, or if binning was specified
                        that.link_needs_update = false;
                        return;
                    }
                    var positions, normals;
                    if (do_clean) {
                        var pn = that.clean_positions_and_normals(bin_size);
                        positions = pn.positions;
                        normals = pn.normals;
                    } else {
                        positions = that.get_positions(geometry.attributes.position.array);
                        normals = that.get_normals(geometry.attributes.normal.array);
                    }
                    geometry.attributes.position.array = positions;
                    geometry.attributes.position.needsUpdate = true;
                    geometry.attributes.normal.array = normals;
                    geometry.attributes.normal.needsUpdate = true;
                    geometry.attributes.color.array = that.get_colors(geometry.attributes.color.array);
                    geometry.attributes.normal.needsUpdate = true;
                    that.link_needs_update = false;
                }
                this.settings.after_run_callback = after_run;
                this.check_update_link = check_update_link;
                return geometry;
            };
            clean_positions_and_normals(normal_binning, truncate) {
                var positions = this.get_positions();
                var normals = this.get_normals();
                var nfloats = positions.length;
                var clean_positions = new Float32Array(nfloats);
                var clean_normals = new Float32Array(nfloats);
                var clean_length = 0;
                var tetrahedron_indices = this.crossing.compact_indices;
                var vertices_per_tetrahedron = this.segments.vertices_per_instance;
                var too_small = this.settings.invalid_coordinate + 1;
                var maxes = null;
                var mins = null;
                for (var i=0; i<tetrahedron_indices.length; i++) {
                    if (tetrahedron_indices[i] < 0) {
                        break;  // sentinel: end of valid tetrahedron indices
                    }
                    var tetrahedron_start = 3 * i * vertices_per_tetrahedron;
                    for (var vj=0; vj<vertices_per_tetrahedron; vj++) {
                        var vertex_start = 3 * vj + tetrahedron_start;
                        if (positions[vertex_start] > too_small) {
                            if (!maxes) {
                                maxes = [];
                                mins = [];
                                for (var k=0; k<3; k++) {
                                    var p = positions[vertex_start + k];
                                    maxes.push(p);
                                    mins.push(p);
                                }
                            }
                            for (var k=0; k<3; k++) {
                                var copy_index = vertex_start + k;
                                var p = positions[copy_index];
                                maxes[k] = Math.max(maxes[k], p);
                                mins[k] = Math.min(mins[k], p)
                                clean_positions[clean_length] = p;
                                clean_normals[clean_length] = normals[copy_index];
                                clean_length++;
                            }
                        }
                    }
                }
                if (normal_binning && (clean_length > 0)) {
                    debugger;
                    // unify geometrically close normal values
                    var key_to_normal = {};
                    var denominators = [];
                    for (var i=0; i<3; i++) {
                        var d = maxes[i] - mins[i];
                        if (d < 1e-17) {
                            d = 1.0
                        }
                        denominators.push(d);
                    }
                    var position_bin_key = function (vertex_index) {
                        var key = 0;
                        var vertex_start = 3 * vertex_index;
                        for (var k=0; k<3; k++) {
                            key = normal_binning * key;
                            var coordinate = clean_positions[vertex_start + k];
                            var k_offset = Math.floor(normal_binning * (coordinate - mins[k])/denominators[k]);
                            key += k_offset;
                        }
                        return key;
                    };
                    var n_vertices = clean_length / 3;
                    var vertex_to_key = {};
                    var key_to_normal_sum = {};
                    for (var vi=0; vi<n_vertices; vi++) {
                        var key = position_bin_key(vi);
                        vertex_to_key[vi] = key;
                        var ns = key_to_normal_sum[key];
                        if (!ns) {
                            ns = [0, 0, 0];
                        }
                        var vertex_start = vi * 3;
                        for (var k=0; k<3; k++) {
                            ns[k] += clean_normals[vertex_start + k];
                        }
                        key_to_normal_sum[key] = ns;
                    }
                    // renormalize
                    for (var k in key_to_normal_sum) {
                        var ns = key_to_normal_sum[k];
                        var n = 0;
                        for (var k=0; k<3; k++) {
                            n += ns[k] * ns[k];
                        }
                        if (n < 1e-10) {
                            n = 1.0;
                        }
                        n = Math.sqrt(n);
                        for (var k=0; k<3; k++) {
                            ns[k] = ns[k] / n;
                        }
                        key_to_normal_sum[k] = ns;
                    }
                    // apply unified normals
                    for (var vi=0; vi<n_vertices; vi++) {
                        var vertex_start = vi * 3;
                        var key = vertex_to_key[vi];
                        var ns = key_to_normal_sum[key];
                        for (var k=0; k<3; k++) {
                            clean_normals[vertex_start + k] = ns[k];
                        }
                    }
                }
                if (truncate) {
                    clean_positions = clean_positions.subarray(0, clean_length);
                    clean_normals = clean_normals.subarray(0, clean_length);
                }
                return {
                    positions: clean_positions,
                    normals: clean_normals,
                    length: clean_length,
                    maxes: maxes,
                    mins: mins,
                }
            }
            set_grid_limits(grid_mins, grid_maxes) {
                this.crossing.set_grid_limits(grid_mins, grid_maxes);
            };
            set_threshold(value) {
                this.settings.threshold = value;
                this.crossing.set_threshold(value);
                // xxxx must be after first run!
                if (this.segments) {
                    this.segments.set_threshold(value);
                }
            };
            get_positions(a) {
                return this.segments.get_positions(a);
            };
            get_normals(a) {
                return this.segments.get_normals(a);
            };
            get_colors(a) {
                return this.segments.get_colors(a);
            };
        };

        return new WebGL2Surfaces3dOpt(options);
    };

    $.fn.webGL2surfaces3d = function (options) {

        // XXXX THIS IS HISTORICAL AND HAS NOT BEEN UPDATED FOR NEW CONVENTIONS XXXX

        // from grid of sample points generate iso-surfacde triangulation.
        class WebGL2Surfaces3d {
            constructor(options) {
                this.settings = $.extend({
                    // default settings:
                    feedbackContext: null,    // the underlying FeedbackContext context to use
                    valuesArray: null,   // the array buffer of values to contour
                    num_rows: null,
                    num_cols: null,
                    num_layers: 1,  // default to "flat"
                    dx: [1, 0, 0],
                    dy: [0, 1, 0],
                    dz: [0, 0, 1],
                    translation: [-1, -1, 0],
                    color: [1, 1, 1],
                    //rasterize: false,
                    threshold: 0,  // value at contour
                    invalid_coordinate: -100000,  // invalidity marker for positions
                    after_run_callback: null,   // call this after each run.
                }, options);
                var s = this.settings;
                this.feedbackContext = s.feedbackContext;
                if (!this.feedbackContext) {
                    throw new Error("Feedback context required.");
                }
                var nvalues = s.valuesArray.length;
                var nvoxels = s.num_rows * s.num_cols * s.num_layers * s.num_blocks;
                if (nvalues != nvoxels) {
                    // for now strict checking
                    throw new Error("voxels " + nvoxels + " don't match values " + nvalues);
                }
                // allocate and load buffer with a fresh name
                this.buffer = this.feedbackContext.buffer()
                this.buffer.initialize_from_array(s.valuesArray);
                var buffername = this.buffer.name;

                this.program = this.feedbackContext.program({
                    vertex_shader: tetrahedra_vertex_shader,
                    fragment_shader: tetrahedra_fragment_shader,
                    feedbacks: {
                        vPosition: {num_components: 3},
                        vNormal: {num_components: 3},
                        vColor: {num_components: 3},
                    },
                })

                // set up input parameters
                var x_offset = 1;
                var y_offset = s.num_cols;
                var z_offset = s.num_cols * s.num_rows;
                var num_instances = nvalues - (x_offset + y_offset + z_offset);

                var inputs = {};
                var add_input = function (ix, iy, iz) {
                    var name = (("a" + ix) + iy) + iz;
                    var dx = [0, x_offset][ix];
                    var dy = [0, y_offset][iy];
                    var dz = [0, z_offset][iz];
                    inputs[name] = {
                        per_vertex: false,
                        num_components: 1,
                        from_buffer: {
                            name: buffername,
                            skip_elements: dx + dy + dz,
                        }
                    }
                };
                add_input(0,0,0);
                add_input(0,0,1);
                add_input(0,1,0);
                add_input(0,1,1);
                add_input(1,0,0);
                add_input(1,0,1);
                add_input(1,1,0);
                add_input(1,1,1);

                const N_TETRAHEDRA = 6;
                const N_TRIANGLES = 2;  
                const N_VERTICES = 3;
                var vertices_per_instance = N_TETRAHEDRA * N_TRIANGLES * N_VERTICES;
                // add vertex count bogus input for Firefox
                var vertexNumArray = new Float32Array(Array.from(Array(vertices_per_instance).keys()));
                this.vertex_num_buffer = this.feedbackContext.buffer()
                this.vertex_num_buffer.initialize_from_array(vertexNumArray);
                inputs["aVertexCount"] = {
                    per_vertex: true,
                    num_components: 1,
                    from_buffer: {
                        name: this.vertex_num_buffer.name,
                    }
                }

                this.runner = this.program.runner({
                    run_type: "TRIANGLES",
                    num_instances: num_instances,
                    vertices_per_instance: vertices_per_instance,
                    rasterize: s.rasterize,
                    uniforms: {
                        uRowSize: {
                            vtype: "1iv",
                            default_value: [s.num_cols],
                        },
                        uColSize: {
                            vtype: "1iv",
                            default_value: [s.num_rows],
                        },
                        uValue: {
                            vtype: "1fv",
                            default_value: [s.threshold],
                        },
                        dx: {
                            vtype: "3fv",
                            default_value: s.dx,
                        },
                        dy: {
                            vtype: "3fv",
                            default_value: s.dy,
                        },
                        dz: {
                            vtype: "3fv",
                            default_value: s.dz,
                        },
                        translation: {
                            vtype: "3fv",
                            default_value: s.translation,
                        },
                        uInvalid: {
                            vtype: "1fv",
                            default_value: [s.invalid_coordinate],
                        },
                    },
                    inputs: inputs,
                });
            };
            run() {
                // may not always need to do re-install uniforms?
                this.runner.install_uniforms();
                this.runner.run();
                var after_run_callback = this.settings.after_run_callback;
                if (after_run_callback) {
                    after_run_callback(this);
                }
            };
            linked_three_geometry (THREE) {
                // create a three.js geometry linked to the current positions feedback array.
                // xxxx only one geometry may be linked at a time.
                // this is a bit convoluted in an attempt to only update attributes when needed.
                var that = this;
                var positions = this.get_positions();
                var normals = this.get_normals();
                var colors = this.get_colors();
                var geometry = new THREE.BufferGeometry();
                geometry.setAttribute( 'position', new THREE.BufferAttribute( positions, 3 ) );
                geometry.setAttribute( 'normal', new THREE.BufferAttribute( normals, 3 ) );
                geometry.setAttribute( 'color', new THREE.BufferAttribute( colors, 3 ) );
                that.link_needs_update = false;
                var after_run = function(that) {
                    debugger;
                    that.link_needs_update = true;
                }
                var check_update_link = function() {
                    // update the geometry positions array in place and mark for update in geometry
                    if (! that.link_needs_update) {
                        // only update upon request and only if needed
                        that.link_needs_update = false;
                        return;
                    }
                    geometry.attributes.position.array = that.get_positions(geometry.attributes.position.array);
                    geometry.attributes.position.needsUpdate = true;
                    geometry.attributes.normal.array = that.get_normals(geometry.attributes.normal.array);
                    geometry.attributes.normal.needsUpdate = true;
                    geometry.attributes.color.array = that.get_colors(geometry.attributes.color.array);
                    geometry.attributes.normal.needsUpdate = true;
                    that.link_needs_update = false;
                }
                this.settings.after_run_callback = after_run;
                this.check_update_link = check_update_link;
                return geometry;
            };
            set_threshold(value) {
                //this.runner.uniforms.uValue.value = [value];
                this.runner.change_uniform("uValue", [value]);
                //this.runner.run();
            };
            get_positions(optionalPreAllocatedArrBuffer) {
                return this.runner.feedback_array(
                    "vPosition",
                    optionalPreAllocatedArrBuffer);
            };
            get_normals(optionalPreAllocatedArrBuffer) {
                return this.runner.feedback_array(
                    "vNormal",
                    optionalPreAllocatedArrBuffer);
            };
            get_colors(optionalPreAllocatedArrBuffer) {
                return this.runner.feedback_array(
                    "vColor",
                    optionalPreAllocatedArrBuffer);
            };
        };
        var tetrahedra_vertex_shader = `#version 300 es
        // triangulate tetrahedral tiling on voxel with vertex values 
        //  a000 .. a111
        // Each voxel is divided into 6 tetrahedra with
        // each tetrahedron split by (up to) 2 triangles.

        // global length of rows, cols
        uniform int uRowSize;
        uniform int uColSize;
        
        // global contour threshold
        uniform float uValue;
        
        // uniform offsets in xyz directions
        uniform vec3 dx, dy, dz, translation;
        
        // invalid value marker
        uniform float uInvalid;
        
        // per mesh function values at voxel corners
        in float a000, a001, a010, a011, a100, a101, a110, a111;

        // Which vertex in which triangle on which tetrahedron?
        //   encodes tetrahedron_number 0..4, triangle_number 0..1, vertex number 0..2
        //   for a total of 5 * 2 * 3 = 30 vertices per "mesh instance".
        //   aVertexCount = tetrahedron_number * 10 + triangle_number * 3 + vertex_number;
        const int N_TETRAHEDRA = 6; // tetrahedra per cube
        const int N_TRIANGLES = 2;  // triangles per tetrahedron
        const int N_VERTICES = 3;   // vertices per triangle
        const int N_CORNERS = 8;    // number of cube corners
        const int N_T_VERTICES = 4; // number of vertices in a tetrahedron

        // bogus vertex attribute required by Firefox (but not Chrome)
        in float aVertexCount;

        // Crossing index is binary integer of form
        //   (triangle_num << 4) || ((fA > v) << 3 || ((fB > v) << 2 || ((fC> v) << 1 || ((fD > v)
        const int N_CROSSING_INDEX = 32;  // 2 ** 5

        // corner offsets
        const vec3 offsets[N_CORNERS] = vec3[] (
            vec3(0.0, 0.0, 0.0),
            vec3(0.0, 0.0, 1.0),
            vec3(0.0, 1.0, 0.0),
            vec3(0.0, 1.0, 1.0),
            vec3(1.0, 0.0, 0.0),
            vec3(1.0, 0.0, 1.0),
            vec3(1.0, 1.0, 0.0),
            vec3(1.0, 1.0, 1.0)
        );

        // vertex indices for tiling tetrahedra per tetrahedron number
        const int A_INDEX = 0;
        const int[N_TETRAHEDRA] B_index = int[] (4, 6, 2, 3, 1, 5) ;
        const int[N_TETRAHEDRA] C_index = int[] (6, 2, 3, 1, 5, 4) ;
        const int D_INDEX = 7;

        // crossing index to triangle vertices endpoint indices
        const int U_left[N_CROSSING_INDEX] = int[] (
            -1, 3, 2, 0, 1, 0, 0, 0, 0, 1, 1, 1, 2, 2, 3,-1,
            -1,-1,-1, 0,-1, 0, 0,-1,-1, 1, 1,-1, 2,-1,-1,-1);
        const int U_right[N_CROSSING_INDEX] = int[] (
            -1, 1, 1, 2, 2, 1, 1, 2, 2, 0, 0, 2, 0, 1, 1,-1,
            -1,-1,-1, 3,-1, 3, 2,-1,-1, 3, 2,-1, 1,-1,-1,-1);
        const int V_left[N_CROSSING_INDEX] = int[] (
            -1, 3, 2, 1, 1, 0, 3, 0, 0, 2, 1, 1, 3, 2, 3,-1,
            -1,-1,-1, 1,-1, 2, 3,-1,-1, 2, 3,-1, 3,-1,-1,-1);
        const int V_right[N_CROSSING_INDEX] = int[] (
            -1, 0, 3, 2, 0, 3, 1, 1, 3, 0, 2, 3, 0, 0, 2,-1,
            -1,-1,-1, 2,-1, 3, 1,-1,-1, 0, 2,-1, 0,-1,-1,-1);
        const int W_left[N_CROSSING_INDEX] = int[] (
            -1, 3, 2, 0, 1, 2, 0, 0, 0, 1, 3, 1, 2, 2, 3,-1,
            -1,-1,-1, 1,-1, 2, 3,-1,-1, 2, 3,-1, 3,-1,-1,-1);
        const int W_right[N_CROSSING_INDEX] = int[] (
            -1, 2, 0, 3, 3, 1, 2, 3, 1, 3, 0, 0, 1, 3, 0,-1,
            -1,-1,-1, 3,-1, 1, 2,-1,-1, 3, 0,-1, 1,-1,-1,-1);

        // feedbacks out
        out vec3 vColor, vPosition, vNormal;

        // debugging
        out float[4] vdump;

        void main() {
            // initially set output point to invalid
            gl_Position = vec4(uInvalid, uInvalid, uInvalid, uInvalid);
            // use the bogus vertexCount parameter so it is not erased by the optimizer
            float grey = aVertexCount / float(N_TETRAHEDRA * N_TRIANGLES * N_VERTICES);
            vColor = vec3(float(gl_VertexID) * 0.01, grey, 0.0);  // temp value for debugging
            vNormal = vec3(0.0, 0.0, 1.0);    // arbitrary initial value

            // size of layer of rows and columns in 3d grid
            int layer_size = uRowSize * uColSize;
            // instance depth of this layer
            int i_depth_num = gl_InstanceID / layer_size;
            // ravelled index in layer
            int i_layer_index = gl_InstanceID - (i_depth_num * layer_size);
            // instance row
            int i_row_num = i_layer_index / uRowSize;
            // instance column
            int i_col_num = i_layer_index - (i_row_num * uRowSize);
            // float versions for calculations
            float layer_num = float(i_depth_num);
            float row_num = float(i_row_num);
            float col_num = float(i_col_num);

            // Dont tile last column which wraps around or last row
            if (i_col_num < (uRowSize - 1) && (i_row_num < (uColSize - 1))) {
                // determine which vertex in which triangle in which tetrahedron to interpolate
                int iVertexCount = gl_VertexID;
                int iTetrahedronNumber = iVertexCount / (N_TRIANGLES * N_VERTICES);
                int iTetIndex = iVertexCount - (N_TRIANGLES * N_VERTICES) * iTetrahedronNumber;
                int iTriangleNumber = iTetIndex / N_VERTICES;
                int iVertexNumber = iTetIndex - (iTriangleNumber * N_VERTICES);
                // offsets of vertices for this tet number
                vec3 t_offsets[N_T_VERTICES] = vec3[](
                    offsets[A_INDEX],
                    offsets[B_index[iTetrahedronNumber]],
                    offsets[C_index[iTetrahedronNumber]],
                    offsets[D_INDEX]
                );
                // weights as array
                float wts[N_CORNERS] = float[](
                    a000, a001, a010, a011, a100, a101, a110, a111);
                // weights of vertices for this tet number
                float t_wts[N_T_VERTICES] = float[](
                    wts[A_INDEX],
                    wts[B_index[iTetrahedronNumber]],
                    wts[C_index[iTetrahedronNumber]],
                    wts[D_INDEX]
                );
                vdump = t_wts;

                // crossing index
                int ci = iTriangleNumber << 1;
                if (t_wts[0] > uValue) { ci = ci + 1; }
                ci = ci << 1;
                if (t_wts[1] > uValue) { ci = ci + 1; }
                ci = ci << 1;
                if (t_wts[2] > uValue) { ci = ci + 1; }
                ci = ci << 1;
                if (t_wts[3] > uValue) { ci = ci + 1; }

                // If U_left[ci] for this corner is negative (invalid index) then there is no such triangle here.
                if (U_left[ci] >= 0) {
                    int SegLs[N_VERTICES] = int[](U_left[ci], V_left[ci], W_left[ci]);
                    int SegRs[N_VERTICES] = int[](U_right[ci], V_right[ci], W_right[ci]);
                    
                    int SegL = SegLs[iVertexNumber];
                    int SegR = SegRs[iVertexNumber];
                    vec3 offsetL = t_offsets[SegL];
                    vec3 offsetR = t_offsets[SegR];
                    
                    float wtL = t_wts[SegL];
                    float wtR = t_wts[SegR];
                    // check denominator is not too small? xxxx
                    float delta = (wtL - uValue) / (wtL - wtR);
                    vec3 combined_offset = ((1.0 - delta) * offsetL) + (delta * offsetR);
                    vec3 vertex = combined_offset + vec3(col_num, row_num, layer_num);
                    vPosition = dx * vertex[0] + dy * vertex[1] + dz * vertex[2] + translation;
                    gl_Position.xyz = vPosition;
                    gl_Position[3] = 1.0;
                    vdump = float[4](vertex[0], vertex[1], vertex[2], delta);

                    // Compute normal for the whole tetrahedron
                    vec3 center = (t_offsets[0] + t_offsets[1] + t_offsets[2] + t_offsets[3])/4.0;
                    vec3 nm = ( 
                        + (t_offsets[0] - center) * (t_wts[0] - uValue) 
                        + (t_offsets[1] - center) * (t_wts[1] - uValue) 
                        + (t_offsets[2] - center) * (t_wts[2] - uValue) 
                        + (t_offsets[3] - center) * (t_wts[3] - uValue) 
                        );
                    float ln = length(nm);
                    if (ln > 1e-12) {
                        vNormal = nm / ln;
                    }
                    vColor = abs(vNormal);  // XXX FOR TESTING ONLY
                }
            }
            vPosition = gl_Position.xyz;
        }
        `;
        
        var tetrahedra_fragment_shader = `#version 300 es
            #ifdef GL_ES
                precision highp float;
            #endif
            in vec3 vColor;
            out vec4 color;
    
            void main() {
                color = vec4(vColor, 1.0);
            }
            `;
        
        return new WebGL2Surfaces3d(options);
    };

    $.fn.webGL2surfaces3d.simple_example = function (container, opt) {
        var gl = $.fn.feedWebGL2.setup_gl_for_example(container);

        var context = container.feedWebGL2({
            gl: gl,
        });
        var valuesArray = new Float32Array([
            0,0,0,
            0,0,0,
            0,0,0,

            0,0,0,
            0,1,0,
            0,0,0,

            0,0,0,
            0,-1,0,
            0,0,0,
        ]);
        var h = 0.5
        var ddz = 0.1
        var init = container.webGL2surfaces3d;
        if (opt) {
            init = container.webGL2surfaces3dopt;
        }
        var contours = init(
            {
                feedbackContext: context,
                valuesArray: valuesArray,
                num_rows: 3,
                num_cols: 3,
                num_layers: 3,
                dx: [h, 0, 0],
                dy: [0, h, 0],
                dz: [ddz, 0.33*ddz, h],
                translation: [-h, -h, -h],
                color: [h, h, h],
                rasterize: true,
                threshold: 0.3,
                // only for "optimized"
                shrink_factor: 0.8,
            }
        );
        // attach an input to change the threshold
        $("<br/>").appendTo(container);
        $("<span>Threshold: </span>").appendTo(container);
        var input = $('<input value="0.3" type="text"></p> ').appendTo(container);
        var dump = $("<div>data dump here</div>").appendTo(container);
        var update = (function () {
            var threshold = + input.val();
            gl.clearColor(0.8, 0.9, 1.0, 1.0);
            gl.clear(gl.COLOR_BUFFER_BIT);
            contours.set_threshold(threshold);
            contours.run();
            var tf = function(x) { return " " + x.toFixed(2)  + " "; };
            var positions = contours.get_positions();
            dump.empty();
            $("<div>" + positions.length + " POSITIONS </div>").appendTo(dump);
            for (var i=0; i<positions.length; i+=3) {
                if (true || positions[i] > -100) {
                    $("<div>" + 
                    tf(positions[i])+ 
                    tf(positions[i+1])+ 
                    tf(positions[i+2])+ 
                    "</div>").appendTo(dump);
                }
            }
            var normals = contours.get_normals();
            $("<br/>").appendTo(dump);
            $("<div>" + normals.length + " NORMALS </div>").appendTo(dump);
            for (var i=0; i<normals.length; i+=3) {
                if (true || normals[i] > -100) {
                    $("<div>" + 
                    tf(normals[i])+ 
                    tf(normals[i+1])+ 
                    tf(normals[i+2])+ 
                    "</div>").appendTo(dump);
                }
            }
        });
        input.change(update);
        update();
        return contours;
    };

})(jQuery);