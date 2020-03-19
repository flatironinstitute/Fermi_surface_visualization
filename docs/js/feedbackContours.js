
/*

JQuery plugin populating 2 and 3 dimensional contours.

Requires nd_frame to be loaded.

Structure follows: https://learn.jquery.com/plugins/basic-plugin-creation/

*/
"use strict";

(function($) {

    $.fn.webGL2crossingPixels = function(options) {
        class WebGL2CrossingPixels {
            constructor(options) {
                this.settings = $.extend({
                    feedbackContext: null,    // the underlying FeedbackContext context to use
                    valuesArray: null,   // the array buffer of values to contour
                    num_rows: null,
                    num_cols: null,
                    num_layers: 1,  // default to "flat"
                    rasterize: false,
                    threshold: 0,  // value at contour
                    // when getting compact arrays
                    // shrink the array sizes by this factor.
                    shrink_factor: 0.2,
                }, options);
                var s = this.settings;
                this.feedbackContext = s.feedbackContext;
                var nvalues = s.valuesArray.length;
                var nvoxels = s.num_rows * s.num_cols * s.num_layers;
                if (nvalues != nvoxels) {
                    // for now strict checking
                    throw new Error("voxels " + nvoxels + " don't match values " + nvalues);
                }
                // allocate and load buffer with a fresh name
                this.buffer = this.feedbackContext.buffer()
                this.buffer.initialize_from_array(s.valuesArray);
                var buffername = this.buffer.name;

                var crossingPixelsVertexShader = `#version 300 es

                // global length of rows
                uniform int uRowSize;

                // global number of columnss
                uniform int uColSize;
                
                // global contour threshold
                uniform float uValue;

                // per array values at pixel corners
                in float aLL, aLR, aUL, aUR;

                // corners feedback
                out vec4 corners;

                // index feedback
                flat out int index;

                void main() {
                    // initially assume the pixed doesn't cross
                    index = -1;
                    corners = vec4(aLL, aLR, aUL, aUR);

                    // size of layer of rows and columns in 3d grid
                    int layer_size = uRowSize * uColSize;
                    // instance depth of this layer
                    int i_depth_num = gl_VertexID / layer_size;
                    // ravelled index in layer
                    int i_layer_index = gl_VertexID - (i_depth_num * layer_size);

                    int i_row_num = i_layer_index/ uRowSize;
                    int i_col_num = i_layer_index - (i_row_num * uRowSize);
                    // Dont tile last column which wraps around rows
                    //vdump = float[](float(uRowSize), col_num, row_num, float(gl_InstanceID));
                    if ((i_col_num < (uRowSize - 1)) && (i_row_num < (uColSize - 1))) {
                        // the pixel doesn't wrap around the edge
                        float m = corners[0];
                        float M = corners[0];
                        for (int i=1; i<4; i++) {
                            float c = corners[i];
                            m = min(c, m);
                            M = max(c, M);
                        }
                        if ((m <= uValue) && (M > uValue)) {
                            // the pixel crosses the threshold
                            index = gl_VertexID;
                        }
                    }
                }
                `;
                this.program = this.feedbackContext.program({
                    vertex_shader: crossingPixelsVertexShader,
                    fragment_shader: this.settings.fragment_shader,
                    feedbacks: {
                        index: {type: "int"},
                        corners: {num_components:4},
                    },
                });

                var x_offset = 1;
                var y_offset = s.num_cols;
                //var z_offset = s.num_cols * s.num_rows;
                var num_pixels = nvalues - (x_offset + y_offset);

                this.runner = this.program.runner({
                    vertices_per_instance: num_pixels,
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
                    },
                    inputs: {
                        aLL: {
                            per_vertex: true,
                            num_components: 1,
                            from_buffer: {
                                name: buffername,
                                skip_elements: 0,
                                element_stride: 1,
                            },
                        },
                        aLR: {
                            per_vertex: true,
                            num_components: 1,
                            from_buffer: {
                                name: buffername,
                                skip_elements: x_offset,
                            },
                        },
                        aUL: {
                            per_vertex: true,
                            num_components: 1,
                            from_buffer: {
                                name: buffername,
                                skip_elements: y_offset,
                            },
                        },
                        aUR: {
                            per_vertex: true,
                            num_components: 1,
                            from_buffer: {
                                name: buffername,
                                skip_elements: x_offset + y_offset,
                            },
                        },
                    }
                });
                this.corners_array = null;
                this.index_array = null;
                this.compact_length = null;
            };
            run() {
                this.runner.install_uniforms();
                this.runner.run();
            };
            get_compacted_feedbacks() {
                this.run();
                var rn = this.runner;
                this.corners_array = rn.feedback_array(
                    "corners",
                    this.corners_array,
                );
                this.index_array = rn.feedback_array(
                    "index",
                    this.index_array,
                );
                if (this.compact_length === null) {
                    // allocate arrays to size limit
                    this.compact_length = Math.floor(
                        this.settings.shrink_factor * this.index_array.length
                    );
                    this.compact_corners = new Float32Array(4 * this.compact_length);
                    this.compact_indices = new Int32Array(this.compact_length);
                }
                // compact the arrays
                this.compact_indices = this.feedbackContext.filter_degenerate_entries(
                    this.index_array, this.index_array, this.compact_indices, 1, -1
                );
                this.compact_corners = this.feedbackContext.filter_degenerate_entries(
                    this.index_array, this.corners_array, this.compact_corners, 4, -1
                );
                return {indices: this.compact_indices, corners: this.compact_corners};
            };
            set_threshold(value) {
                this.runner.change_uniform("uValue", [value]);
            };
        };
        return new WebGL2CrossingPixels(options);
    };

    $.fn.webGL2crossingPixels.example = function (container) {
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
        var crossing = container.webGL2crossingPixels({
            feedbackContext: context,
            valuesArray: valuesArray,
            num_rows: 3,
            num_cols: 3,
            num_layers: 3,  // default to "flat"
            threshold: 0.5,
            shrink_factor: 0.3,
        });
        var compacted = crossing.get_compacted_feedbacks();
        var indices = compacted.indices;
        var corners = compacted.corners;
        var ci = 0;
        for (var i=0; i<indices.length; i++) {
            $("<br/>").appendTo(container);
            $("<span> " + indices[i] + " </span>").appendTo(container);
            for (var j=0; j<4; j++) {
                $("<span> " + corners[ci] + " </span>").appendTo(container);
                ci ++;
            }
        }
    };

    $.fn.webGL2SegmentPixels = function (options) {
        class WebGL2SegmentPixels {
            constructor(options) { 
                this.settings = $.extend({
                    feedbackContext: null,
                    // array of indices (from crossing pixels)
                    indices: null,
                    // array of corners (from crossing pixels)
                    corners: null,
                    num_rows: null,
                    num_cols: null,
                    dx: [1, 0, 0],
                    dy: [0, 1, 0],
                    dz: [0, 0, 1],
                    translation: [-1, -1, 0],
                    color: [1, 1, 1],
                    rasterize: false,
                    threshold: 0,  // value at contour
                    invalid_coordinate: -100000,  // invalidity marker for positions
                }, options);
                var s = this.settings;
                this.feedbackContext = s.feedbackContext;
                
                // allocate and load buffers with a fresh name
                this.index_buffer = this.feedbackContext.buffer()
                this.index_buffer.initialize_from_array(s.indices);
                this.corner_buffer = this.feedbackContext.buffer()
                this.corner_buffer.initialize_from_array(s.corners);
                // add vertex count bogus input for Firefox
                var vertexNumArray = new Float32Array([0,1,2,3])
                this.vertex_num_buffer = this.feedbackContext.buffer()
                this.vertex_num_buffer.initialize_from_array(vertexNumArray);

                this.program = this.feedbackContext.program({
                    vertex_shader: segments_vertex_shader,
                    feedbacks: {
                        vPosition: {num_components: 3},
                    },
                });

                this.runner = this.program.runner({
                    run_type: "LINES",
                    num_instances: s.indices.length,
                    vertices_per_instance: 4,
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
                    inputs: {
                        index: {
                            per_vertex: false,
                            num_components: 1,
                            type: "int",
                            from_buffer: {
                                name: this.index_buffer.name,
                            },
                        },
                        corners: {
                            per_vertex: false,
                            num_components: 4,
                            from_buffer: {
                                name: this.corner_buffer.name,
                            },
                        },
                        aVertexCount: {   // bogus attribute required by Firefox
                            per_vertex: true,
                            num_components: 1,
                            from_buffer: {
                                name: this.vertex_num_buffer.name,
                            }
                        }
                    }
                });
            };
            run() {
                this.runner.install_uniforms();
                this.runner.run();
            };
            set_threshold(value) {
                this.runner.change_uniform("uValue", [value]);
            };
            get_positions(optionalPreAllocatedArrBuffer) {
                return this.runner.feedback_array(
                    "vPosition",
                    optionalPreAllocatedArrBuffer);
            };
        };
        var segments_vertex_shader = `#version 300 es

        // global length of rows
        uniform int uRowSize;

        // global number of columnss
        uniform int uColSize;
        
        // global contour threshold
        uniform float uValue;

        // xxxxx add divisor for multiple contours...
        
        // uniform offsets in x,y,z directions, translation, line color
        uniform vec3 dx, dy, dz, translation, color;

        // invalid value marker
        uniform float uInvalid;

        // per mesh pixel index
        in int index;

        // per mesh corner values
        in vec4 corners;

        // per vertex -- which vertex? 0,1 on first triangle or 2,3 on second
        in float aVertexCount;  // bogus attribute needed by Firefox

        // feedback output position of vertex (or degenerate)
        out vec3 vPosition;
        // square corner offsets
        const vec2 offsets[4] = vec2[](
            vec2(0.0,0.0),
            vec2(1.0,0.0),
            vec2(0.0,1.0),
            vec2(1.0,1.0)
        );

        // crossing index to segment endpoint indices
        //                            000 001 010 011 100 101 110 111
        const int Seg1Left[8] = int[]( -1,  1,  0,  2,  2,  0,  1, -1);
        const int Seg1Right[8]= int[]( -1,  2,  1,  0,  0,  1,  2, -1);
        const int Seg2Left[8] = int[]( -1,  2,  1,  0,  0,  1,  2, -1);
        const int Seg2Right[8]= int[]( -1,  0,  2,  1,  1,  2,  0, -1);

        void main() {
            // fake use of aVertexCount to prevent optimizer from eliminating it.
            vPosition[2] = aVertexCount; 
            gl_Position = vec4(uInvalid, uInvalid, uInvalid, uInvalid);
            vPosition = gl_Position.xyz;

            // size of layer of rows and columns in 3d grid
            int layer_size = uRowSize * uColSize;
            // instance depth of this layer
            int i_depth_num = index / layer_size;
            // ravelled index in layer
            int i_layer_index = index - (i_depth_num * layer_size);

            int i_row_num = i_layer_index/ uRowSize;
            int i_col_num = i_layer_index - (i_row_num * uRowSize);

            // Dont tile last column which wraps around rows
            if ((index>=0) && (i_col_num < (uRowSize - 1)) && (i_row_num < (uColSize - 1))) {
                float row_num = float(i_row_num);
                float col_num = float(i_col_num);
                float depth_num = float(i_depth_num);
                // determine which vertex in which triangle to interpolate
                int iVertexCount = gl_VertexID;
                int iTriangleNumber = iVertexCount / 2;
                int iVertexNumber = iVertexCount - (iTriangleNumber * 2);
                int iT1 = iTriangleNumber + 1;
                vec2 triangle_offsets[3] = vec2[](
                    offsets[0],
                    offsets[iT1],
                    offsets[3]
                );
                float triangle_wts[3] = float[](
                    corners[0],
                    corners[iT1],
                    corners[3]
                );
                // crossing index
                int ci = 0;
                for (int i=0; i<3; i++) {
                    ci = ci << 1;
                    if (triangle_wts[i] > uValue) {
                        ci += 1;
                    }
                }
                if (Seg1Left[ci] >= 0) {
                    int SegLs[2] = int[](Seg1Left[ci], Seg2Left[ci]);
                    int SegRs[2] = int[](Seg1Right[ci], Seg2Right[ci]);
                    int SegL = SegLs[iVertexNumber];
                    int SegR = SegRs[iVertexNumber];
                    vec2 offsetL = triangle_offsets[SegL];
                    vec2 offsetR = triangle_offsets[SegR];
                    float wtL = triangle_wts[SegL];
                    float wtR = triangle_wts[SegR];
                    // check denominator is not too small? xxxx
                    float delta = (wtL - uValue) / (wtL - wtR);
                    vec2 combined_offset = ((1.0 - delta) * offsetL) + (delta * offsetR);
                    //combined_offset = offsetL;
                    vec2 vertex = combined_offset + vec2(col_num, row_num);
                    vPosition = 
                        dx * vertex[0] + 
                        dy * vertex[1] + 
                        dz * depth_num + 
                        translation;
                    gl_Position.xyz = vPosition;
                    //vColor = abs(normalize(vec3(vertex)));  // XXX FOR TESTING ONLY
                    gl_Position[3] = 1.0;
                }
            }
            //vPosition = gl_Position.xyz;
        }
        `;
        return new WebGL2SegmentPixels(options);
    };

    $.fn.webGL2SegmentPixels.example = function (container) {
        var gl = $.fn.feedWebGL2.setup_gl_for_example(container);

        var context = container.feedWebGL2({
            gl: gl,
        });
        var cornersArray = new Float32Array([
            0, 0, 0, 1,
            0, 0, 1, 0,
            0, 1, 0, 0,
            1, 0, 0, 0,
            -1, -1, -1, -1, -1,
            -1, -1, -1, -1, -1,
        ]);
        var indexArray = new Int32Array([
            9, 10, 12, 13, -1, -1,
        ]);
        var segments = container.webGL2SegmentPixels({
            feedbackContext: context,
            indices: indexArray,
            corners: cornersArray,
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

    $.fn.webGL2contours2dopt = function (options) {
        // "optimized contours" by truncating buffer sizes
        // which may result in some data omission in dense cases.
        class WebGL2Contour2d {
            constructor(options) {
                this.settings = $.extend({
                    // default settings:
                    shrink_factor: 0.1, // how much to shrink buffers
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
                    rasterize: false,
                    threshold: 0,  // value at contour
                    invalid_coordinate: -100000,  // invalidity marker for positions
                    after_run_callback: null,   // call this after each run.
                }, options);
                var s = this.settings;
                this.feedbackContext = s.feedbackContext;
                var container = $(this.feedbackContext.canvas);
                if (!this.feedbackContext) {
                    throw new Error("Feedback context required.");
                }
                var nvalues = s.valuesArray.length;
                var nvoxels = s.num_rows * s.num_cols * s.num_layers;
                if (nvalues != nvoxels) {
                    // for now strict checking
                    throw new Error("voxels " + nvoxels + " don't match values " + nvalues);
                }
                this.crossing = container.webGL2crossingPixels({
                    feedbackContext: this.feedbackContext,
                    valuesArray: s.valuesArray,
                    num_rows: s.num_rows,
                    num_cols: s.num_cols,
                    num_layers: s.num_layers,  // default to "flat"
                    threshold: s.threshold,
                    shrink_factor: s.shrink_factor,
                    // never rasterize the crossing pixels
                });
                // initialize segmenter upon first run.
                this.segments = null; 
            };
            run () {
                var s = this.settings;
                var compacted = this.crossing.get_compacted_feedbacks();
                if (!this.segments) {
                    var container = $(this.feedbackContext.canvas);
                    this.segments = container.webGL2SegmentPixels({
                        feedbackContext: this.feedbackContext,
                        indices: compacted.indices,
                        corners: compacted.corners,
                        num_rows: s.num_rows,
                        num_cols: s.num_cols,
                        rasterize: s.rasterize,
                        dx: s.dx,
                        dy: s.dy,
                        dz: s.dz,
                        translation: s.translation,
                        threshold: s.threshold,
                    });
                } else {
                    // reset buffer content
                    this.segments.index_buffer.copy_from_array(
                        compacted.indices
                    );
                    this.segments.corner_buffer.copy_from_array(
                        compacted.corners
                    );
                }
                this.segments.run();
                var after_run_callback = this.settings.after_run_callback;
                if (after_run_callback) {
                    after_run_callback(this);
                }
                //var positions = segments.get_positions();
            };
            linked_three_geometry (THREE) {
                // create a three.js geometry linked to the current positions feedback array.
                // xxxx only one geometry may be linked at a time.
                var that = this;
                var positions = this.get_positions();
                var geometry = new THREE.BufferGeometry();
                geometry.setAttribute( 'position', new THREE.BufferAttribute( positions, 3 ) );
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
                    that.link_needs_update = false;
                }
                this.settings.after_run_callback = after_run;
                this.check_update_link = check_update_link;
                return geometry;
            };
            set_threshold(value) {
                debugger;
                this.crossing.set_threshold(value);
                // xxxx must be after first run!
                if (this.segments) {
                    this.segments.set_threshold(value);
                }
            };
            get_positions() {
                return this.segments.get_positions();
            }
        };

        return new WebGL2Contour2d(options);
    };

    $.fn.webGL2contours2d = function (options) {
        // from grid of sample points generate contour line segments

        class WebGL2Contour2d {
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
                    rasterize: false,
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
                var nvoxels = s.num_rows * s.num_cols * s.num_layers;
                if (nvalues != nvoxels) {
                    // for now strict checking
                    throw new Error("voxels " + nvoxels + " don't match values " + nvalues);
                }
                // allocate and load buffer with a fresh name
                this.buffer = this.feedbackContext.buffer()
                this.buffer.initialize_from_array(s.valuesArray);
                var buffername = this.buffer.name;

                this.program = this.feedbackContext.program({
                    vertex_shader: contour_vertex_shader,
                    fragment_shader: contour_fragment_shader,
                    feedbacks: {
                        vPosition: {num_components: 3},
                    },
                });
                var x_offset = 1;
                var y_offset = s.num_cols;
                //var z_offset = s.num_cols * s.num_rows;
                var num_instances = nvalues - (x_offset + y_offset);
                var vertices_per_instance = 4; // 2 endpoints each for 2 triangles.
                // add vertex count bogus input for Firefox
                var vertexNumArray = new Float32Array(Array.from(Array(vertices_per_instance).keys()));
                this.vertex_num_buffer = this.feedbackContext.buffer()
                this.vertex_num_buffer.initialize_from_array(vertexNumArray);

                this.runner = this.program.runner({
                    run_type: "LINES",
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
                    inputs: {
                        aLL: {
                            per_vertex: false,
                            num_components: 1,
                            from_buffer: {
                                name: buffername,
                                skip_elements: 0,
                                element_stride: 1,
                            },
                        },
                        aLR: {
                            per_vertex: false,
                            num_components: 1,
                            from_buffer: {
                                name: buffername,
                                skip_elements: x_offset,
                            },
                        },
                        aUL: {
                            per_vertex: false,
                            num_components: 1,
                            from_buffer: {
                                name: buffername,
                                skip_elements: y_offset,
                            },
                        },
                        aUR: {
                            per_vertex: false,
                            num_components: 1,
                            from_buffer: {
                                name: buffername,
                                skip_elements: x_offset + y_offset,
                            },
                        },
                        aVertexCount: {   // bogus attribute required by Firefox
                            per_vertex: true,
                            num_components: 1,
                            from_buffer: {
                                name: this.vertex_num_buffer.name,
                            }
                        }
                    }
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
                var that = this;
                var positions = this.get_positions();
                var geometry = new THREE.BufferGeometry();
                geometry.setAttribute( 'position', new THREE.BufferAttribute( positions, 3 ) );
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
                    that.link_needs_update = false;
                }
                this.settings.after_run_callback = after_run;
                this.check_update_link = check_update_link;
                return geometry;
            };
            set_threshold(value) {
                this.runner.change_uniform("uValue", [value]);
            };
            get_positions(optionalPreAllocatedArrBuffer) {
                return this.runner.feedback_array(
                    "vPosition",
                    optionalPreAllocatedArrBuffer);
            };
        };
        var contour_vertex_shader = `#version 300 es
            // triangulate contour segments on pixel with vertex values aLL, aLR, aUL, aUR
            //
            // UL   V0   UR
            //  .---*---.
            //  | T0   /|
            //  |     / |
            //  * V1 /  |
            //  |   /   * V3
            //  |  /    |
            //  | /  T1 |
            //  .---*---.
            // LL   V2   LR
            //
            // output vertices are segment (V1, V0) crossing triangle T0, if exists
            //    and segment (V2, V3) crossing triangle T1 if exists.
            // The segments can be oriented opposite any of the 3 triangle vertices.

            // global length of rows
            uniform int uRowSize;

            // global number of columnss
            uniform int uColSize;
            
            // global contour threshold
            uniform float uValue;

            // xxxxx add divisor for multiple contours...
            
            // uniform offsets in x,y,z directions, translation, line color
            uniform vec3 dx, dy, dz, translation, color;
            
            // invalid value marker
            uniform float uInvalid;
            
            // per mesh function values at pixel corners
            in float aLL, aLR, aUL, aUR;

            // per vertex -- which vertex? 0,1 on first triangle or 2,3 on second
            in float aVertexCount;  // bogus attribute needed by Firefox

            // feedbacks out
            out vec3 vColor, vPosition;

            // debugging
            out float[4] vdump;

            // square corner offsets
            const vec2 offsets[4] = vec2[](
                vec2(0.0,0.0),
                vec2(1.0,0.0),
                vec2(0.0,1.0),
                vec2(1.0,1.0)
            );

            // crossing index to segment endpoint indices
            //                            000 001 010 011 100 101 110 111
            const int Seg1Left[8] = int[]( -1,  1,  0,  2,  2,  0,  1, -1);
            const int Seg1Right[8]= int[]( -1,  2,  1,  0,  0,  1,  2, -1);
            const int Seg2Left[8] = int[]( -1,  2,  1,  0,  0,  1,  2, -1);
            const int Seg2Right[8]= int[]( -1,  0,  2,  1,  1,  2,  0, -1);

            void main() {
                // initially set output point to invalid
                gl_Position = vec4(uInvalid, uInvalid, uInvalid, uInvalid);
                vPosition = gl_Position.xyz;
                vColor = vec3(aVertexCount * 0.1, 0.0, 0.0);

                // size of layer of rows and columns in 3d grid
                int layer_size = uRowSize * uColSize;
                // instance depth of this layer
                int i_depth_num = gl_InstanceID / layer_size;
                // ravelled index in layer
                int i_layer_index = gl_InstanceID - (i_depth_num * layer_size);

                int i_row_num = i_layer_index/ uRowSize;
                int i_col_num = i_layer_index - (i_row_num * uRowSize);
                // Dont tile last column which wraps around rows
                //vdump = float[](float(uRowSize), col_num, row_num, float(gl_InstanceID));
                if ((i_col_num < (uRowSize - 1)) && (i_row_num < (uColSize - 1))) {
                    float row_num = float(i_row_num);
                    float col_num = float(i_col_num);
                    float depth_num = float(i_depth_num);
                    vdump = float[](aLL, aLR, aUL, aUR);
                    // determine which vertex in which triangle to interpolate
                    int iVertexCount = gl_VertexID;
                    int iTriangleNumber = iVertexCount / 2;
                    int iVertexNumber = iVertexCount - (iTriangleNumber * 2);
                    int iT1 = iTriangleNumber + 1;
                    vec2 triangle_offsets[3] = vec2[](
                        offsets[0],
                        offsets[iT1],
                        offsets[3]
                    );
                    float wts[4] = float[](aLL, aLR, aUL, aUR);
                    float triangle_wts[3] = float[](
                        wts[0],
                        wts[iT1],
                        wts[3]
                    );
                    //vdump = wts;
                    // crossing index
                    int ci = 0;
                    for (int i=0; i<3; i++) {
                        ci = ci << 1;
                        if (triangle_wts[i] > uValue) {
                            ci += 1;
                        }
                    }
                    if (Seg1Left[ci] >= 0) {
                        int SegLs[2] = int[](Seg1Left[ci], Seg2Left[ci]);
                        int SegRs[2] = int[](Seg1Right[ci], Seg2Right[ci]);
                        int SegL = SegLs[iVertexNumber];
                        int SegR = SegRs[iVertexNumber];
                        vec2 offsetL = triangle_offsets[SegL];
                        vec2 offsetR = triangle_offsets[SegR];
                        float wtL = triangle_wts[SegL];
                        float wtR = triangle_wts[SegR];
                        // check denominator is not too small? xxxx
                        float delta = (wtL - uValue) / (wtL - wtR);
                        vec2 combined_offset = ((1.0 - delta) * offsetL) + (delta * offsetR);
                        //combined_offset = offsetL;
                        vec2 vertex = combined_offset + vec2(col_num, row_num);
                        vPosition = 
                            dx * vertex[0] + 
                            dy * vertex[1] + 
                            dz * depth_num + 
                            translation;
                        gl_Position.xyz = vPosition;
                        //vColor = abs(normalize(vec3(vertex)));  // XXX FOR TESTING ONLY
                        gl_Position[3] = 1.0;
                    }
                }
                vPosition = gl_Position.xyz;
            }
            `;
        
        var contour_fragment_shader = `#version 300 es
            #ifdef GL_ES
                precision highp float;
            #endif
            in vec3 vColor;
            out vec4 color;
    
            void main() {
                color = vec4(vColor, 1.0);
            }
            `;
        
        return new WebGL2Contour2d(options);
    };

    $.fn.webGL2contours2d.simple_example = function (container, opt) {
        var gl = $.fn.feedWebGL2.setup_gl_for_example(container);

        var context = container.feedWebGL2({
            gl: gl,
        });
        var valuesArray = new Float32Array([
            0,0,0,
            0,1,0,
            0,0,0,

            0,0,0,
            0,2,0,
            0,0,0,

            0,0,0,
            0,0,0,
            0,0,0,
        ]);
        var h = 0.5;
        var ddz = 0.1;
        var init = container.webGL2contours2d;
        if (opt) {
            init = container.webGL2contours2dopt;
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
                shrink_factor: 0.5,
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
            for (var i=0; i<positions.length; i+=4) {
                if (true || (positions[i] > -100)) {
                    $("<div>" + 
                    tf(positions[i])+ 
                    tf(positions[i+1])+ 
                    tf(positions[i+2])+ 
                    tf(positions[i+3])+ "</div>").appendTo(dump);
                }
            }
            if (opt) {
                var compacted = contours.crossing.get_compacted_feedbacks();
                var indices = compacted.indices;
                var corners = compacted.corners;
                var ci = 0;
                for (var ii=0; ii<indices.length; ii++) {
                    $("<br/>").appendTo(dump);
                    $("<span> " + indices[ii] + " </span>").appendTo(dump);
                    for (var j=0; j<4; j++) {
                        $("<span> " + corners[ci] + " </span>").appendTo(dump);
                        ci ++;
                    }
                }
            }
        });
        input.change(update);
        update();
        return contours;
    };

})(jQuery);