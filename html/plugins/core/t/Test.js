dojo.provide("plugins.core.Test");

dojo.require("plugins.core.Common");

dojo.declare("plugins.core.Test",
    plugins.core.Common,
//[ plugins.core.Common, plugins.core.form.Template ],
    {
        //Path to the template of this widget. 
        cssFiles : [ "file1.css" ],

        constructor : function ()
        {

            var div = document.createElement('div');
            document.body.appendChild(div);
            div.innerHTML = "<h1>Loaded</h1>";
        }

        //,


        //loadCSS : function ()
        //{
        //    
        //    
        //}

    }
);

