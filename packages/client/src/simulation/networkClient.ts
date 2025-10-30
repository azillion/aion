export class NetworkClient {
  private ws: WebSocket | null = null;
  public onData: ((data: ArrayBuffer) => void) | null = null;
  public onConnect: (() => void) | null = null;

  public connect(url: string): Promise<void> {
    return new Promise((resolve, reject) => {
      this.ws = new WebSocket(url);
      this.ws.binaryType = 'arraybuffer';

      this.ws.onopen = () => {
        console.log(`[NetworkClient] Connected to ${url}`);
        this.onConnect?.();
        resolve();
      };

      this.ws.onmessage = (event: MessageEvent) => {
        if (event.data instanceof ArrayBuffer) {
          this.onData?.(event.data);
          return;
        }
        // Some browsers deliver Blob even with binaryType set; normalize to ArrayBuffer
        if (event.data instanceof Blob) {
          event.data.arrayBuffer().then((buf) => {
            this.onData?.(buf);
          }).catch((err) => {
            console.error('[NetworkClient] Failed to read Blob data:', err);
          });
          return;
        }
        // Fallback: try to handle typed arrays
        if (ArrayBuffer.isView(event.data)) {
          const view = event.data as ArrayBufferView;
          this.onData?.(view.buffer);
          return;
        }
        console.warn('[NetworkClient] Unsupported message type', typeof event.data);
      };

      this.ws.onerror = (event) => {
        console.error('[NetworkClient] WebSocket error:', event);
        reject(new Error('WebSocket connection failed.'));
      };

      this.ws.onclose = () => {
        console.log('[NetworkClient] Connection closed.');
      };
    });
  }

  public disconnect(): void {
    this.ws?.close();
  }
}


