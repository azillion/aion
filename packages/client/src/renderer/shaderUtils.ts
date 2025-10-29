import type { GPUDevice } from '@webgpu/types';

const includeRegex = /#include\s+"([^"]+)"/g;

export async function createShaderModule(
  device: GPUDevice,
  label: string,
  code: string,
  shaderIncludes: Record<string, string>
): Promise<GPUShaderModule> {
  
  const processedCode = code.replace(includeRegex, (match, path) => {
    const includeContent = shaderIncludes[path];
    if (includeContent === undefined) {
      console.error(`Shader include not found: ${path}`);
      return `// ERROR: SHADER INCLUDE "${path}" NOT FOUND`;
    }
    return includeContent;
  });

  return device.createShaderModule({
    label,
    code: processedCode,
  });
}


