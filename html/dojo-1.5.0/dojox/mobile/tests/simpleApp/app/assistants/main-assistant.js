dojo.provide("MainAssistant");
dojo.require("dojox.mobile.app.SceneAssistant");

dojo.declare("MainAssistant", dojox.mobile.app.SceneAssistant, {

  setup: function(){

    var appInfoNode = this.controller.query(".appInfoArea")[0];

    appInfoNode.innerHTML = 
      "This app has the following info: \n" 
        + dojo.toJson(dojox.mobile.app.info, true)

  },

  activate: function(){


  }

});