// WHAT WE'RE GETTING
dojo.provide("plugins.admin.ConfirmDialogue");

// NEED THESE FOR THE WIDGETS IN THE TEMPLATE
dojo.require("dojox.layout.FloatingPane");
dojo.require("dijit.TitlePane");
dojo.require("dijit.layout.TabContainer");
dojo.require("dijit.form.ComboBox"); 
dojo.require("dijit.form.Button"); 
dojo.require("dojo.parser"); // scan page for widgets

// INHERITED CLASS
dojo.require("plugins.core.form.Template");

// HERE IT IS
dojo.declare(
    "plugins.admin.ConfirmDialogue",
	plugins.core.form.Template,
	{
        //Path to the template of this widget. 
        templatePath: dojo.moduleUrl("plugins", "admin/templates/confirmDialogue.html"),

        // Calls dijit._Templated.widgetsInTemplate
        widgetsInTemplate : true,

        // PARENT WIDGET
        parentWidget    : null,

        // VALUES
        title           :   "Please confirm",
        message         :   "Do you really want to delete this application?",
        option1         :   "Cancel",
        option2         :   "Delete it",

        widgets : [  ],

        // Any initialization code would go here in the constructor.
		// plugins.report.Template and its superclasses dijit._Widget and
        // dijit._Templated do not have parameters in their constructors, so
        // there wouldn't be any multiple-inheritance complications
        // if you were to include some paramters here.
        constructor : function (args)
        {

            // SET this.args
            this.title = args.title;
            this.message = args.message;
            this.option1 = args.option1;
            this.option2 = args.option2;

            // SET PARENT AND GRANDPARENT IDS
            this.parentWidget = args.parentWidget;

        },


        postCreate : function (args)
        {




        }
    }

);
