
// http://localhost:8080/misc/volume/nih/surfaces/example.html

class SurfaceController {
    constructor(json_data, target_id, value_id, column_id, row_id, layer_id) {
        var that = this;
        this.json_data = json_data;
        var s = json_data;
        this.$container = $("#" + target_id);
        this.$value = $("#" + value_id);
        this.$columns = $("#" + column_id);
        this.$rows = $("#" + row_id);
        this.$layers = $("#" + layer_id);
        this.container = this.$container[0];
        this.$container.html("JSON loaded.");
        this.context = this.$container.feedWebGL2({});
        this.valuesArray = new Float32Array(s.values);
        var middle = 0.5 * (s.min_value + s.max_value);
        this.grid_mins = [-1, -1, -1];
        this.grid_maxes = [s.num_columns+1, s.num_rows+1, s.num_layers+1];
        this.surfaces = this.$container.webGL2surfaces3dopt(
              {
                  feedbackContext: this.context,
                  valuesArray: this.valuesArray,
                  num_rows: s.num_rows,
                  num_cols: s.num_columns,
                  num_layers: s.num_layers,
                  rasterize: false,
                  threshold: middle,
                  shrink_factor: 0.1,  // how much to shrink the arrays
                  grid_min: this.grid_mins,
                  grid_max: this.grid_maxes,
              }
        );
        this.$value.empty();
        this.slider = $('<div></div>').appendTo(this.$value);
        this.slider.width(this.$value.width() * 0.8)
        this.slider_readout = $('<div>readout</div>').appendTo(this.$value)

        /*
        this.$layers.empty();
        this.layers_slider = $('<div></div>').appendTo(this.$layers);
        this.layers_slider.width(this.$layers.width() * 0.8)
        */

        var update = function () {
            var threshold = + that.slider.slider("option", "value");
            that.slider_readout.html("value: " + threshold.toFixed(4))
            that.surfaces.set_threshold(threshold);

            /*
            var layer_limits = that.layers_slider.slider("option", "values");
            that.grid_mins[2] = layer_limits[0];
            that.grid_maxes[2] = layer_limits[1];
            */
            for (var name in that.grid_sliders) {
                var slider = that.grid_sliders[name];
                var index = slider.grid_index;
                var limits = slider.slider("option", "values");
                that.grid_mins[index] = limits[0];
                that.grid_maxes[index] = limits[1];
            }
            that.surfaces.set_grid_limits(that.grid_mins, that.grid_maxes);

            that.surfaces.run();
            
            if (that.scene_initialized) {
                that.update_scene();
            } else {
                that.initialize_scene();
            }
        };

        this.grid_sliders = {
            columns: this.make_grid_slider({
                name: "columns",
                index: 0,
                $target: this.$columns,
                update: update,
                max: s.num_columns,
            }),
            rowss: this.make_grid_slider({
                name: "rows",
                index: 1,
                $target: this.$rows,
                update: update,
                max: s.num_rows,
            }),
            layers: this.make_grid_slider({
                name: "layers",
                index: 2,
                $target: this.$layers,
                update: update,
                max: s.num_layers,
            }),
        }

        /*
        this.layers_slider.slider({
            range: true,
            min: -1,
            max: s.num_layers+1,
            values: [0, s.num_layers],
            slide: update,
            change: update,
        });
        */

        this.slider.slider({
            min: s.min_value,
            max: s.max_value,
            value: middle,
            step: 0.01 * (s.max_value - s.min_value),
            slide: update,
            change: update,
        });
        this.scene_initialized = false;

        update();
    };
    make_grid_slider(options) {
        options.$target.html(options.name);
        var slider = $("<div>  </div>").appendTo(options.$target);
        slider.width(options.$target.width() * 0.8)
        slider.slider({
            range: true,
            min: -1,
            max: options.max + 1,
            values: [0, options.max],
            slide: options.update,
            change: options.update,
        });
        slider.grid_index = options.index;
        return slider;
    };
    initialize_scene() {
        var that = this;
        var container = this.container;
        var $container = this.$container;
        var canvas = document.createElement( 'canvas' ); 
        var context = canvas.getContext( 'webgl2', { alpha: false } ); 
        this.renderer = new THREE.WebGLRenderer( { canvas: canvas, context: context } );
        //renderer = new THREE.WebGLRenderer();
        this.renderer.setPixelRatio( window.devicePixelRatio );
        this.renderer.setSize( $container.width() * 0.99, $container.height()*0.80);
        //renderer.setSize( window.innerWidth, window.innerHeight );
        this.renderer.outputEncoding = THREE.sRGBEncoding;
        $container.empty();
        container.appendChild( this.renderer.domElement );
        this.scene = new THREE.Scene();
        this.camera = new THREE.PerspectiveCamera( 45, $container.width()/$container.height(), 1, 10000 );
        this.camera.position.set( 0, 0, -3 );
        this.geometry = this.surfaces.linked_three_geometry(THREE);
        var material = new THREE.MeshNormalMaterial( {  } );
        material.side = THREE.DoubleSide;
        this.mesh = new THREE.Mesh( this.geometry,  material );
        this.scene.add( this.mesh );
        this.orbitControls = new THREE.OrbitControls(this.camera, this.renderer.domElement);
        this.orbitControls.userZoom = false;
        this.clock = new THREE.Clock();
        this.update_scene();
        this.scene_initialized = true;
        var animate = function () {
            that.surfaces.check_update_link();
            var delta = that.clock.getDelta();
            that.orbitControls.update(delta);
            that.renderer.render( that.scene, that.camera );
            requestAnimationFrame( animate );
        };
        animate();
    };
    update_scene() {
      // update is automagic because there is an animation loop
    };
};
