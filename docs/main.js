
// globals for easy debugging
var json_data, name_to_controller, wire_frame;
var normals, velocities, collector;

class BufferCollector {
    constructor(level_controllers, runner) {
        this.runner = runner;
        this.level_controllers = level_controllers;
        // assuming level controllers have been run at least once.
        var vertex_buffer_size = 0;
        for (var name in level_controllers) {
            var c = level_controllers[name];
            c.buffer_start_index = vertex_buffer_size;
            vertex_buffer_size += c.vertices.length;
        }
        this.vertex_buffer_size = vertex_buffer_size;
        this.vertices = new Float32Array(vertex_buffer_size);
        this.normals = new Float32Array(vertex_buffer_size);
        this.velocity_colors = new Float32Array(vertex_buffer_size);
        this.wire_colors = new Float32Array(vertex_buffer_size);
        // set up wire colors only once
        this.copy_arrays("wire_colors");
        this.copy_all_mutable_arrays();
        // set up geometries
        var wf = new THREE.BufferGeometry();
        this.wireframe_geometry = wf;
        wf.setAttribute( 'position', new THREE.BufferAttribute( this.vertices, 3 ) );
        wf.setAttribute( 'color', new THREE.BufferAttribute( this.wire_colors, 3 ) );
        var n = new THREE.BufferGeometry();
        this.normal_geometry = n;
        n.setAttribute( 'position', new THREE.BufferAttribute( this.vertices, 3 ) );
        n.setAttribute( 'normal', new THREE.BufferAttribute( this.normals, 3 ) );
        var v = new THREE.BufferGeometry();
        this.velocity_geometry = v;
        v.setAttribute( 'position', new THREE.BufferAttribute( this.vertices, 3 ) );
        v.setAttribute( 'color', new THREE.BufferAttribute( this.velocity_colors, 3 ) );
    };
    recalibrate()  {
        this.runner();
        this.copy_all_mutable_arrays();
        this.update_geometries();
    }
    update_geometries() {
        var wfa = this.wireframe_geometry.attributes;
        wfa.position.array = this.vertices;
        wfa.position.needsUpdate = true;
        // don't update wirecolor
        var na = this.normal_geometry.attributes;
        na.position.array = this.vertices;
        na.position.needsUpdate = true;
        na.normal.array = this.normals;
        na.normal.needsUpdate = true;
        var va = this.velocity_geometry.attributes;
        va.position.array = this.vertices;
        va.position.needsUpdate = true;
        va.color.array = this.velocity_colors;
        va.color.needsUpdate = true;
    }
    copy_arrays(array_name) {
        var level_controllers = this.level_controllers;
        var target_array = this[array_name];
        for (var name in level_controllers) {
            var c = level_controllers[name];
            var start = c.buffer_start_index;
            var source_array = c[array_name];
            target_array.set(source_array, start)
        }
    };
    copy_all_mutable_arrays() {
        this.copy_arrays("vertices");
        this.copy_arrays("normals");
        this.copy_arrays("velocity_colors");
    };
};

var default_layer_colors = [
    [1.0, 0.5, 0.5],
    [0.0, 0.8, 0.8],
    [7.0, 0.7, 0.3],
    [0.7, 0.2, 0.7],
];

class EnergyLevelController {
    // aggregate controller for normal and velocity for an energy level (0 or 1)
    constructor(level, local_data, json_data) {
        var valuesArray = new Float32Array(local_data.values);
        this.colors = new Float32Array(local_data.colors);
        this.wire_color = local_data.wire_color 
            || default_layer_colors[level]
            || default_layer_colors[0];
        // unattached dom container (?)
        this.container = $("<div/>");
        this.context = this.container.feedWebGL2({});
        var middle = 0.5 * (json_data.min_value + json_data.max_value);
        var surface_options = {
            feedbackContext: this.context,
            valuesArray: valuesArray,
            num_rows: json_data.num_rows,
            num_cols: json_data.num_columns,
            num_layers: json_data.num_layers,
            rasterize: false,
            threshold: middle,
            shrink_factor: 0.4,
        }
        this.surfaces = this.container.webGL2surfaces3dopt(surface_options);
        this.wire_colors = null;
        this.vertices = null;
        this.normals = null;
        this.velocity_colors = null;
        this.buffer_start_index = null;
    };
    run() {
        var surfaces = this.surfaces;
        surfaces.set_threshold(json_data.threshold);
        surfaces.set_grid_limits(json_data.grid_mins, json_data.grid_maxes);
        surfaces.run();
        this.vertices = surfaces.get_positions(this.vertices);
        this.normals = surfaces.get_normals(this.normals);
        if (!this.wire_colors) {
            // initialize colors arrays
            var ln = this.vertices.length;
            var w = new Float32Array(ln);
            var c = this.wire_color;
            for (var i=0; i<ln; i++) {
                w[i] = c[i % 3];
            }
            this.wire_colors = w;
            this.velocity_colors = new Float32Array(ln);
        }
        // XXXXX TEMPORARY DEBUG VALUE FOR VELOCITY_COLORS
        //this.velocity_colors = this.wire_colors;
        this.velocity_colors = surfaces.colorization(this.colors, this.velocity_colors);
    };
};

class NormalsController {
    // controls the displays colored by triangle normals
    constructor(identifier, collector) {
        this.collector = collector;
        this.$container = $("#" + identifier);
        this.container = this.$container[0];
        this.canvas = document.createElement( 'canvas' ); 
        this.context = this.canvas.getContext( 'webgl2', { alpha: false } ); 
        this.renderer = new THREE.WebGLRenderer( { canvas: this.canvas, context: this.context } );
        this.renderer.setPixelRatio( window.devicePixelRatio );
        this.renderer.setSize( this.$container.width() * 0.99, this.$container.height()*0.95);
        //renderer.setSize( window.innerWidth, window.innerHeight );
        this.renderer.outputEncoding = THREE.sRGBEncoding;
        this.container.appendChild( this.renderer.domElement );

        this.scene = new THREE.Scene();
        this.scene.add( new THREE.AmbientLight( 0x444444 ) );
        this.camera = new THREE.PerspectiveCamera( 45, this.$container.width()/this.$container.height(), 1, 10000 );
        this.camera.position.set( 0, 0, 4 );
        this.geometry = this.get_geometry();
        this.material = this.get_material();
        this.material.side = THREE.DoubleSide;
        this.mesh = new THREE.Mesh( this.geometry,  this.material );
        this.scene.add( this.mesh );
        this.renderer.render( this.scene, this.camera );
    };
    get_geometry() {
        return this.collector.normal_geometry;
    };
    get_material() {
        return new THREE.MeshNormalMaterial( {  } );
    };
    sync_camera(other_camera) {
        var camera = this.camera;
        camera.position.x = other_camera.position.x;
        camera.position.y = other_camera.position.y;
        camera.position.z = other_camera.position.z;
        camera.lookAt(new THREE.Vector3(0, 0, 0));
        this.renderer.render( this.scene, this.camera );
    };
};

class VelocityController extends NormalsController {
    // controls displays colorized by velocity
    get_geometry() {
        return this.collector.velocity_geometry;
    };
    get_material() {
        return new THREE.MeshBasicMaterial( {vertexColors: THREE.VertexColors, side: THREE.DoubleSide} );
    };
};

class WireFrameController extends NormalsController {
    get_geometry() {
        return this.collector.wireframe_geometry;
    };
    get_material() {
        var material = new THREE.MeshBasicMaterial({ vertexColors: THREE.VertexColors });
        material.wireframe = true;
        return material;
    };
    orbit(normals, velocities) {
        var that = this;
        this.orbitControls = new THREE.OrbitControls(this.camera, this.renderer.domElement);
        this.orbitControls.userZoom = false;
        this.clock = new THREE.Clock();

        var animate = function () {
            //that.surfaces.check_update_link();
            var delta = that.clock.getDelta();
            that.orbitControls.update(delta);
            that.renderer.render( that.scene, that.camera );
            normals.sync_camera(that.camera);
            velocities.sync_camera(that.camera);
            requestAnimationFrame( animate );
        };
        animate();
    };
};

var set_up = function(data) {
    json_data = data;  // store in global
    // set up the sliders
    var value_slider = $("#value_slider");
    var value_readout = $("#value_readout")
    var m = json_data.min_value;
    var M = json_data.max_value;
    var update_value = function () {
        var v = value_slider.slider("option", "value");
        value_readout.html("" + v.toFixed(4));
        json_data.threshold = v;
        if (collector) {
            collector.recalibrate();
        }
    };
    value_slider.empty();
    value_slider.slider({
        min: m,
        max: M,
        step: 0.01 * (M - m),
        value: 0.5 * (M + m),
        slide: update_value,
        change: update_value,
    });
    update_value();
    json_data.grid_mins = [0, 0, 0];
    json_data.grid_maxes = [json_data.num_columns, json_data.num_rows, json_data.num_layers];
    var col_slider = set_up_dim_slider("col_slider", json_data.num_columns, 0);
    var row_slider = set_up_dim_slider("row_slider", json_data.num_rows, 1);
    var layer_slider = set_up_dim_slider("layer_slider", json_data.num_layers, 2);
    // set up the level controllers
    var level_controllers = {};
    for (var level=0; level<json_data.num_levels; level++) {
        var level_name = "E" + level;
        var level_data = json_data[level_name];
        level_controllers[level_name] = new EnergyLevelController(level, level_data, json_data);
    }
    var run_levels = function () {
        for (var name in level_controllers) {
            level_controllers[name].run();
        }
    };
    run_levels();
    collector = new BufferCollector(level_controllers, run_levels);
    normals = new NormalsController("normals", collector);
    velocities = new VelocityController("velocities", collector);
    wireframe = new WireFrameController("wireframe", collector);
    wireframe.orbit(normals, velocities);
};

set_up_dim_slider = function(container, dim, index) {
    var $container = $("#"+container);
    $container.empty();
    $("<div>" + container + "</div>").appendTo($container);
    var slider = $("<div></div>").appendTo($container);
    var step = Math.max(0.01 * dim, 1);
    var update = function () {
        var limits = slider.slider("option", "values");
        json_data.grid_mins[index] = limits[0];
        json_data.grid_maxes[index] = limits[1];
        if (collector) {
            collector.recalibrate();
        }
    };
    slider.slider({
        range: true,
        min: -1,
        max: dim+1,
        step: step,
        values: [0, dim],
        slide: update,
        change: update,
    });
    return slider;
};

var load_file = function(json_file_url) {
    $.getJSON(json_file_url, set_up).fail(on_load_failure);
};

var on_load_failure = function() {
    alert("Could not load local JSON data.\n" +
            "You may need to run a web server to avoid cross origin restrictions.")
};
