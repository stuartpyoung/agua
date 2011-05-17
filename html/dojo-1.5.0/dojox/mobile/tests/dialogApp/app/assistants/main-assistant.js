dojo.provide("MainAssistant");
dojo.require("dojox.mobile.app.SceneAssistant");

dojo.declare("MainAssistant", dojox.mobile.app.SceneAssistant, {

  setup: function(){

    // Instantiate widgets in the template HTML.
    this.controller.parse();

    var appInfoNode = this.controller.query(".appInfoArea")[0];

    appInfoNode.innerHTML = 
      "This app has the following info: \n" 
        + dojo.toJson(dojox.mobile.app.info, true);

    function handleChoose(value){
      appInfoNode.innerHTML = "Value selected: " + value;

    }


    var controller = this.controller; 


    dojo.connect(dijit.byId("btn1"), "onClick", function(){
      controller.showAlertDialog({
        title: "First Dialog",
        text: "This is a simple text message",
        onChoose: handleChoose,

        buttons: [
          {
            label: "Tap Me!",
            value: "tapped",
            "class": "mblBlueButton"
          }
        ]
      })
    });

    dojo.connect(dijit.byId("btn2"), "onClick", function(){
      controller.showAlertDialog({
        title: "Second Dialog",
        text: "These two buttons return different values, 'value one' and 'value two'",
        onChoose: handleChoose,

        buttons: [
          {
            label: "Im Am Button 1",
            value: "value one",
            "class": "mblBlueButton"
          },
          {
            label: "Im Am Button 2",
            value: "value two",
            "class": "mblBlueButton"
          }
        ]
      })
    });

  },

  activate: function(){


  }

});