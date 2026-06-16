## Pixi and the VS Code Python debugger
Automatic Pixi shell activation can interfere with the Python debugger in VS Code.

As per [this solution](https://github.com/renan-r-santos/pixi-code/issues/21#issuecomment-3413521907), prevent auto-activation of virtual environments in new shells.

```JSON
\\ settings.json
{
  "python.terminal.activateEnvironment": false
}
```
