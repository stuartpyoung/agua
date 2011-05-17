dojo.provide("plugins.form.ValidationTextarea");
dojo.require("dijit.form.SimpleTextarea");
dojo.require("dijit.form.ValidationTextBox");

dojo.declare(
    "plugins.form.ValidationTextarea",
    [dijit.form.ValidationTextBox,dijit.form.SimpleTextarea],
    {
        invalidMessage: "This field is required",

        regExp: "(.|\\s)*",

        preamble : function(args)
        {
            //
            this.invalidMessage = args.invalidMessage;
        },


        constructor : function(args)
        {
            //
            this.invalidMessage = args.invalidMessage;
            this.promptMessage = args.promptMessage;
        },

        postCreate: function() {

            // SAVE INVALID MESSAGE
            var tempInvalidMessage = this.invalidMessage;
            var tempPromptMessage = this.promptMessage;

            this.inherited(arguments);

            // RESTORE INVALID MESSAGEA
            this.invalidMessage = tempInvalidMessage;
            this.promptMessage = tempPromptMessage;
        },

        validate: function() {

            this.inherited(arguments);

            if (arguments.length==0) this.validate(false);
        },

        onFocus: function() {
            if (!this.isValid()) {
                this.displayMessage(this.getErrorMessage());
            }
        },

        onBlur: function() {
            this.validate(false);
        }
     }
);
