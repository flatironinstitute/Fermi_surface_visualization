
// Require jQuery only if needed.
if (!global.jQuery) {
  global.jQuery = require('jquery');
}

// The plugins install themselves into the global jQuery object
require("./feedWebGL");
require("./feedbackContours");
require("./feedbackSurfaces");
require("./feedbackMatrix");

function feedWebGL_is_loaded() {
  return true;
}

export default feedWebGL_is_loaded;
