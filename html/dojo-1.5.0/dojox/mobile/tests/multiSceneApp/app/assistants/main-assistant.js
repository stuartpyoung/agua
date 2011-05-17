dojo.provide("MainAssistant");
dojo.require("dojox.mobile.app.SceneAssistant");

dojo.declare("MainAssistant", dojox.mobile.app.SceneAssistant, {

  setup: function(){

    var appInfoNode = this.controller.query(".appInfoArea")[0];

    appInfoNode.innerHTML = 
      "This app has the following info: \n" 
        + dojo.toJson(dojox.mobile.app.info, true).split("\t").join("  ");

    // Instantiate widgets in the template HTML.
    this.controller.parse();

    var _this = this;
    var launcher = dijit.byId("secondSceneLauncher");

    this.connect(launcher, "onClick", function(){

      _this.controller.stageController.pushScene("second", "Came from Main Scene");
    });
  },

  activate: function(data){

    var node = this.controller.query(".inputData")[0];
    if(data) {
      node.innerHTML = "Scene got the data: " + data;
    } else {
      node.innerHTML = "Scene did not receive data";
    }

  }

});